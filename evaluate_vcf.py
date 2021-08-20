from pysam import VariantFile

import sys
import time
import vcf2hgvs.vcfconverter as vcfconverter
import orjson

if(len(sys.argv) != 3):
    print(f'USAGE: {sys.argv[0]} <VCF_PATH> <JSON_PATH>')
    sys.exit(1)

bcf_in = VariantFile(sys.argv[1])
json_out = open(sys.argv[2], "wb")

start = time.perf_counter()
print('- parsing VCF and creating HGVS data:')

for rec in bcf_in.fetch():
    
    position = vcfconverter.vcf_to_position(rec)    
    json_out.write(orjson.dumps(position, option=orjson.OPT_APPEND_NEWLINE))

print(f'- finished in {time.perf_counter() - start} seconds')
json_out.close()