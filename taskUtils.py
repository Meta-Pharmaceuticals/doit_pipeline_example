import GEOparse
import wget
import tarfile
from pathlib import Path

def download_geo_sup_file(gsenum, outdir):
    outdir.mkdir(parents=True, exist_ok=True)

    gse = GEOparse.get_GEO(geo=gsenum, destdir=outdir)
    for gsm_name, gsm in gse.gsms.items():
        print("Name: ", gsm_name)
        print("Metadata:",)
        for key, value in gsm.metadata.items():
            if key == 'supplementary_file_1':
                urls = value
        for url in urls:
            if url != None:
                tarFileName = wget.download(url, out = str(outdir))
                file = tarfile.open(outdir / tarFileName)
                file.extractall(outdir)
                file.close()