"""Check project-relative static assets, JavaScript syntax and sharing hygiene."""
from pathlib import Path
from html.parser import HTMLParser
import re,subprocess,urllib.parse
ROOT=Path(__file__).resolve().parents[1];WEB=ROOT/'web'

class Links(HTMLParser):
    def handle_starttag(self,tag,attrs):
        for key,value in attrs:
            if key not in ['src','href'] or not value:continue
            parsed=urllib.parse.urlparse(value)
            if parsed.scheme or value.startswith('#'):continue
            assert not value.startswith('/'),f'Root-relative URL breaks project Pages: {value}'
            path=WEB/urllib.parse.unquote(parsed.path)
            assert path.exists(),f'Missing local web asset: {value}'

Links().feed((WEB/'index.html').read_text())
for path in WEB.rglob('*'):
    if path.suffix in ['.js','.mjs']:
        subprocess.run(['node','--check',str(path)],check=True)
        text=path.read_text()
        for dependency in re.findall(r'''(?:from\s*|import\s*)['"](\.[^'"]+)['"]''',text):
            assert (path.parent/dependency).exists(),f'Missing module: {path.name} -> {dependency}'
for directory in [WEB,ROOT/'toolkit',ROOT/'scripts']:
    for path in directory.rglob('*'):
        if not path.is_file() or path.suffix in ['.zip','.npz','.png','.pyc'] or '__pycache__' in path.parts:continue
        if path.resolve()==Path(__file__).resolve():continue
        text=path.read_text(errors='replace')
        # Internal Site identity and private photo access links do not belong in the public package.
        assert 'appgprj_' not in text and 'oai-library://' not in text and 'photos.google.com/share/' not in text, f'Private service reference in {path}'
assert (WEB/'vendor/THREE-LICENSE.txt').exists()
assert not (ROOT/'.openai').exists()
print('PASS: static assets, project-relative paths, JavaScript syntax and distribution checks.')
