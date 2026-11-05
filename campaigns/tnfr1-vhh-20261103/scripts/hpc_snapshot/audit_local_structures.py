"""Inspect deposited numbering and ligand contacts locally, without inference."""
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
AA = dict(zip('ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL'.split(), 'ARNDCQEGHILKMFPSTWYV'))
sequence = json.loads((ROOT / 'sources/P19438.json').read_text())['sequence']['value']

def atoms(pdb, chain):
    rows = []
    for line in (ROOT / 'sources' / pdb).read_text().splitlines():
        if line.startswith('ATOM  ') and line[21] == chain and line[16] in (' ', 'A'):
            if line[76:78].strip() == 'H':
                continue
            rows.append((int(line[22:26]), AA[line[17:20]], line[12:16].strip(),
                         tuple(float(line[s:s+8]) for s in (30, 38, 46))))
    return rows

reference = atoms('9S3I.pdb', 'A')
ca = [a for a in reference if a[2] == 'CA']
for n, aa, _, _ in ca:
    assert aa == sequence[n+29-1], (n, aa)
receptor, ligand = atoms('1TNR.pdb', 'R'), atoms('1TNR.pdb', 'A')
contacts = {}
for n, aa, _, xyz in receptor:
    assert aa == sequence[n+29-1], (n, aa)
    d2 = min(sum((a-b)**2 for a, b in zip(xyz, q[3])) for q in ligand)
    if d2 <= 4.5**2:
        contacts[n+29] = min(contacts.get(n+29, 100), d2**0.5)
result = {
    '9S3I_A': {'deposited_sequence': 'GP + UniProt40-201',
               'author_to_uniprot_offset': 29, 'observed_CA_count': len(ca),
               'observed_uniprot_positions': [r[0]+29 for r in ca],
               'all_observed_residues_match_current_canonical': True},
    '1TNR': {'ligand': 'lymphotoxin-alpha / TNF-beta, not TNF-alpha',
             'method': 'minimum heavy-atom distance <=4.5 A within deposited R:A pair; not full symmetry-expanded biological interface',
             'contacts_uniprot': [{'position': n, 'aa': sequence[n-1], 'distance_A': round(d, 3)}
                                  for n, d in sorted(contacts.items())],
             'interpretation': 'Structural contact candidates, not experimentally established energetic hotspots'},
}
(ROOT / 'reports/local_structure_audit.json').write_text(json.dumps(result, indent=2))
print(json.dumps(result, indent=2))
