"""Verify the compact package against its recorded source hashes."""
from pathlib import Path
import hashlib
import json

ROOT = Path(__file__).resolve().parents[1]

def audit_sources(root=ROOT):
    provenance = root/'provenance'
    manifest = json.loads((provenance/'source_manifest.json').read_text(encoding='utf-8'))
    for name, expected in manifest['files'].items():
        path = root/name
        assert path.is_file(), f'Missing package source: {name}'
        actual = hashlib.sha256(path.read_bytes()).hexdigest()
        assert actual == expected, f'Source changed since packaging: {name}'
    result = {'count': len(manifest['files']), 'status': 'passed',
              'layout': 'single shared function library'}
    (provenance/'source_snapshot_integrity.json').write_text(
        json.dumps(result, indent=2), encoding='utf-8')
    return result

if __name__ == '__main__':
    print(json.dumps(audit_sources(), indent=2))
