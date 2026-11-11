"""Build the standalone toolkit download using the Python standard library only."""
from pathlib import Path
import argparse,zipfile
ROOT=Path(__file__).resolve().parents[1]

def build(output=None):
    output=Path(output) if output else ROOT/'web/Burial_GPR_Toolkit.zip'
    output.parent.mkdir(parents=True,exist_ok=True)
    with zipfile.ZipFile(output,'w',zipfile.ZIP_DEFLATED) as archive:
        for path in sorted((ROOT/'toolkit').rglob('*')):
            if path.is_file() and '__pycache__' not in path.parts and path.suffix!='.pyc':
                archive.write(path,'Burial_GPR_Toolkit/'+str(path.relative_to(ROOT/'toolkit')))
        for name in ['LICENSE','THIRD_PARTY_NOTICES.md','CITATION.cff']:
            archive.write(ROOT/name,'Burial_GPR_Toolkit/'+name)
        for path in sorted((ROOT/'web').rglob('*')):
            if not path.is_file() or path.suffix=='.zip':continue
            target='Burial_GPR_Toolkit/web/'+str(path.relative_to(ROOT/'web'))
            if path.name=='index.html':
                text=path.read_text().replace('href="Burial_GPR_Toolkit.zip" download','href="../README.md"').replace('Download Python toolkit ↓','Read local instructions')
                archive.writestr(target,text)
            else:archive.write(path,target)
    with zipfile.ZipFile(output) as archive:
        bad=archive.testzip()
        if bad:raise RuntimeError('Corrupt ZIP member: '+bad)
    print(f'Built {output.relative_to(ROOT) if output.is_relative_to(ROOT) else output} ({output.stat().st_size:,} bytes)')
    return output

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path)
    build(parser.parse_args().output)
