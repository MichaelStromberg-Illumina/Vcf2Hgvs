from dataclasses import dataclass

import sys
import hgvs.parser
import hgvs.validator
import hgvs.exceptions
import hgvs.dataproviders.uta
import hgvs.assemblymapper
import hgvs.normalizer
import vcf2hgvs.accessions as accessions

hp = hgvs.parser.Parser()
hdp = hgvs.dataproviders.uta.connect()
vr = hgvs.validator.Validator(hdp=hdp)
hn = hgvs.normalizer.Normalizer(hdp)
am = hgvs.assemblymapper.AssemblyMapper(hdp, assembly_name='GRCh37', alt_aln_method='splign', normalize=True, replace_reference=True)

@dataclass
class Transcript(object):
    hgvsc: str
    hgvsp: str

@dataclass
class Variant(object):
    vid: str
    hgvsg: str
    transcripts: list[Transcript]

@dataclass
class Position(object):
    variants: list[Variant]

def hgvsg_to_transcripts(var_g):
    
    transcripts = am.relevant_transcripts(var_g)
    hgvs_transcripts = []
    
    for transcript in transcripts:
        
        if('/' in transcript):
            print(f'WARNING: Skipping {transcript}.')
            continue
            
        try:
            var_t = am.g_to_t(var_g, transcript)
        except hgvs.exceptions.HGVSInvalidIntervalError as e:
            print(f'WARNING: {e}')
            continue
        except hgvs.exceptions.HGVSUnsupportedOperationError as e:
            print(f'ERROR: {e}')
            continue
        except hgvs.exceptions.HGVSDataNotAvailableError as e:
            message = str(e)
            
            # skip transcripts where the alignments are incomplete
            if(message.startswith('Alignment is incomplete')):
                print(f'WARNING: {e}')
                continue;
            
            # skip transcripts that have warnings about non-adjacent exons
            if('are not adjacent' in message):
                print(f'WARNING: {e}')
                continue
            
            print(f'ERROR: {e}')
            sys.exit(1)
        
        try:
            var_p = am.t_to_p(var_t)
        except NotImplementedError as e:
            print(f'WARNING: {e}')
            var_p = None
        except IndexError as e:
            print(f'WARNING: {e}')
            var_p = None
        
        hgvs_transcript = Transcript(str(var_t), str(var_p))
        hgvs_transcripts.append(hgvs_transcript)
    
    return hgvs_transcripts

def vcf_to_position(rec):
    
    variants = []
    
    # convert the chromosome to an accession
    accession = accessions.to_accession(rec.chrom, "GRCh37")
    
    # evaluate each alt allele
    for alt in rec.alts:
        
        # create the VID
        vid = '%s-%d-%s-%s' % (rec.chrom, rec.pos, rec.ref, alt)
        
        # everything can be abstracted as a delins
        hgvs_g = '%s:g.%s_%sdel%sins%s' % (accession, str(rec.pos), str(rec.pos + (len(rec.ref) - 1)), rec.ref, alt)
        var_g = hp.parse_hgvs_variant(hgvs_g)
        
        # check if this validates properly
        try:
            vr.validate(var_g)
        except hgvs.exceptions.HGVSError as e:
            print(f'ERROR: {e}')
            sys.exit()
        
        # normalize the HGVS g. notation
        var_g_norm = hn.normalize(var_g)
        
        # evaluate each transcript
        transcripts = hgvsg_to_transcripts(var_g_norm)
        
        variant = Variant(vid, str(var_g_norm), transcripts)
        variants.append(variant)
    
    position = Position(variants)
    
    return position