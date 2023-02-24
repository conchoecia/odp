# fasta-parser
A lightweight Python module to read and write sequence records in [FASTA format](https://en.wikipedia.org/wiki/FASTA_format).

The design was inspired by the utility of BioPython’s SeqIO, which supports many sequence formats. This repo focuses only on FASTA records. It is faster than BioPython and can handle compressed FASTA files (gzip, bzip2, zip, Zstandard).

## Requirements
Python >= 3.8

### Optional dependencies
This module has no Python package dependencies (i.e., it only uses the standard library). However, if you work with Zstandard-compressed FASTA files, you will need the following package:

* zstandard  >= 0.18.0
    > Install with: `pip install zstandard`

## Quick Start
Typical usage is to read a FASTA file and loop over the record(s).

```python
import fasta

for record in fasta.parse('test/test.fasta'):
    print(record.id, len(record), record.seq[:10], record.desc)
```

Output:

```
NP_002433.1  362   METDAPQPGL   RNA-binding protein Musashi homolog 1 [Homo sapiens]
ENO94161.1    79   MKLLISGLGP   RRM domain-containing RNA-binding protein
sequence     292   MKLSKIALMM
```

## Usage
This module contains the `Record` class representing a FASTA sequence record and the `parse()` function to read FASTA records from a file.

### Record object
Record is an object that contains information on a FASTA sequence record, including id, description, and the sequence itself.

```python
import fasta

record = fasta.Record(
    id='NP_950171.2', 
    seq='MEEEAETEEQQRFSYQQRLKAAVHYTVGCLCEEVALDKEMQFSKQTIAAISELTFRQCENFAKDLEMFASICRKRQE',
    desc='APITD1-CORT protein isoform 2 [Homo sapiens]'
)

print(record.id)            # NP_950171.2
print(record.desc)          # APITD1-CORT protein isoform 2 [Homo sapiens]
print(record.seq)           # MEEEAE..
print(record.description)   # >NP_950171.2 G APITD1-CORT protein isoform 2 [Homo sapiens]
print(len(record))          # 77
print('EEEA' in record)     # True
```

By default, the sequence line in the FASTA format is wrapped to 70 characters. You can provide the line length. Use zero (or None) for no wrapping.

```python
print(record)
# >NP_950171.2 APITD1-CORT protein isoform 2 [Homo sapiens]
# MEEEAETEEQQRFSYQQRLKAAVHYTVGCLCEEVALDKEMQFSKQTIAAISELTFRQCENFAKDLEMFAS
# ICRKRQE

print(record.format(wrap=30))
# >NP_001382951.1 G protein subunit gamma 5 [Homo sapiens]
# MEEEAETEEQQRFSYQQRLKAAVHYTVGCL
# CEEVALDKEMQFSKQTIAAISELTFRQCEN
# FAKDLEMFASICRKRQE

print(record.format(wrap=None))
# >NP_950171.2 APITD1-CORT protein isoform 2 [Homo sapiens]
# MEEEAETEEQQRFSYQQRLKAAVHYTVGCLCEEVALDKEMQFSKQTIAAISELTFRQCENFAKDLEMFASICRKRQE
```

### parse
The `parse()` function is a generator to read FASTA records as `Record` objects one by one from a file (plain FASTA or compressed using gzip or bzip2). Because only one record is created at a time, very little memory is required.

```python
import fasta

for record in fasta.parse('test/test.fasta.gz'):
    print(record.id)
```

For some tasks you may need to have a reusable access to the records. For this purpose, you can use the built-in Python `list()` function to turn the iterator into a list:

```python
import fasta

records = list(fasta.parse('test/test.fasta.gz'))
print(records[0].id)   # First record
print(records[-1].id)  # Last record
```

Another common task is to index your records by sequence identifier. Use `to_dict()` to turn a Record iterator (or list) into a dictionary.

```python
import fasta

records = fasta.to_dict(fasta.parse('test/test.fasta.gz'))
print(records['NP_002433.1'])   # Use any record id
```

## Test
You can run tests to ensure that the module works as expected.

```
./test.py
```

## License

[GNU General Public License, version 3](https://www.gnu.org/licenses/gpl-3.0.html)