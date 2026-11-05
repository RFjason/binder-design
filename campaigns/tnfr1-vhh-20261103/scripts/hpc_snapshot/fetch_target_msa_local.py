"""One ColabFold MSA request/status check per call; sends the public target only."""
import hashlib
import io
import json
from pathlib import Path
import sys
import tarfile
import time

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT.parent / 'IL7RA_20261102' / '.vendor'))
import requests

HOST = 'https://api.colabfold.com'
STATE = ROOT / 'reports/target_msa_ticket.json'
TARGET = json.loads((ROOT / 'sources/P19438.json').read_text())['sequence']['value'][:211]
HEADERS = {'User-Agent': 'TNFR1-Changjiang-Research/1.0'}

def records(text):
    header, parts = None, []
    for line in text.replace('\x00', '').splitlines():
        if line.startswith('>'):
            if header is not None: yield header, ''.join(parts)
            header, parts = line[1:], []
        elif header is not None:
            parts.append(line.strip())
    if header is not None: yield header, ''.join(parts)

if STATE.exists():
    state = json.loads(STATE.read_text())
    assert state['target_sha256'] == hashlib.sha256(TARGET.encode()).hexdigest()
    if state.get('a3m_sha256'):
        print(json.dumps(state))
        raise SystemExit(0)
    response = requests.get(HOST + '/ticket/' + state['id'], headers=HEADERS, timeout=30)
    response.raise_for_status()
    ticket = response.json()
    state['status'] = ticket['status']
else:
    response = requests.post(HOST + '/ticket/msa', data={'q': '>101\n' + TARGET + '\n', 'mode': 'env'}, headers=HEADERS, timeout=30)
    response.raise_for_status()
    ticket = response.json()
    state = {**ticket, 'host': HOST, 'mode': 'env', 'templates': False, 'pairing': False,
             'target_length': len(TARGET), 'binder_sent': False,
             'target_sha256': hashlib.sha256(TARGET.encode()).hexdigest(),
             'created_utc': time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime()),
             'protocol_reference': 'https://github.com/sokrypton/ColabFold/blob/main/colabfold/colabfold.py'}
    if 'id' not in state:
        raise RuntimeError(json.dumps(state))

if state['status'] == 'COMPLETE':
    result = requests.get(HOST + '/result/download/' + state['id'], headers=HEADERS, timeout=40)
    result.raise_for_status()
    archive_path = ROOT / 'sources/TNFR1_target_mmseqs2_env.tar.gz'
    archive_path.write_bytes(result.content)
    entries = []
    source_members = []
    with tarfile.open(fileobj=io.BytesIO(result.content), mode='r:gz') as archive:
        members = {m.name: m for m in archive.getmembers()}
        for name in ('uniref.a3m', 'bfd.mgnify30.metaeuk30.smag30.a3m'):
            match = next((m for key, m in members.items() if key.lstrip('./') == name), None)
            if match is not None:
                source_members.append(name)
                entries.extend(records(archive.extractfile(match).read().decode()))
    assert entries and entries[0][1] == TARGET, 'First MSA row must equal the competition target exactly'
    for header, sequence in entries:
        aligned = ''.join(c for c in sequence if not c.islower() and c != '.')
        assert len(aligned) == len(TARGET), (header, len(aligned))
    folder = ROOT / 'inputs'
    folder.mkdir(exist_ok=True)
    path = folder / 'TNFR1_competition_M1_T211.a3m'
    path.write_text(''.join('>' + h + '\n' + s + '\n' for h, s in entries))
    state.update(a3m_path=str(path), a3m_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
                 rows=len(entries), unique_aligned_rows=len({''.join(c for c in s if not c.islower() and c != '.') for _, s in entries}),
                 result_archive_sha256=hashlib.sha256(result.content).hexdigest(), source_members=source_members)
STATE.write_text(json.dumps(state, indent=2))
print(json.dumps(state), flush=True)
