"""Reading and writing FASTA files.

Copyright (c) 2022 by 
Andrzej Zielezinski (a.zielezinski@gmail.com)
Maciej Michalczyk (mccv99@gmail.com)

https://github.com/aziele/fasta-parser
"""

import bz2      # standard python library
import gzip     # standard python library
import io       # standard python library
import pathlib  # standard python library
import typing   # standard python library
import zipfile  # standard python library

# Zstandard
try:
    import zstandard
    HAS_ZSTANDARD = True
except ImportError:
    HAS_ZSTANDARD = False


class Record:
    """Object representing a FASTA (aka Pearson) record.

    Attributes:
        id  (str)         : Sequence identifier
        seq (str)         : Sequence
        description (str) : Description line (defline)
    """

    def __init__(self, id: str, seq: str, desc: typing.Optional[str] = None):
        """Creates a Record.

        Example:
            >>> record = Record(id='NP_055309.2', 
            ...                 seq='MRELEAKAT',
            ...                 desc='TNRC6A')
            >>> print(record)
            >NP_055309.2 TNRC6A
            MRELEAKAT
        """
        self.id = id
        self.seq = seq
        self.desc = desc

    @property
    def description(self) -> str:
        """Returns a description line (defline) of FASTA record.

        Example:
            >>> record = Record(id='NP_055309.2', 
            ...                 seq='MRELEAKAT',
            ...                 desc='TNRC6A')
            >>> print(record.description)
            >NP_055309.2 TNRC6A
        """
        lst = [f'>{self.id}']
        if self.desc:
            lst.append(f'{self.desc}')
        return " ".join(lst)

    def __iter__(self):
        """Iterates over the letters in the sequence.

        Example:
            >>> record = Record(id='NP_055309.2', 
            ...                 seq='MRELEAKAT',
            ...                 desc='TNRC6A')
            >>> for amino_acid in record:
            ...     print(amino_acid)
            M
            R
            E
            L
            E

            This is equivalent to iterating over the sequence directly:
            >>> for amino_acid in record.seq:
            ...     print(amino_acid)
            M
            R
            E
            L
            E
        """
        return iter(self.seq)

    def __contains__(self, char):
        """Implements the 'in' keyword to search the sequence.

        Example:
            >>> record = Record(id='NP_055309.2', 
            ...                 seq='MRELEAKAT',
            ...                 desc='TNRC6A')
            >>> print('M' in record)
            True
        """
        return char in self.seq

    def __str__(self):
        """Returns the record as a string in the FASTA format.

        Example:
            >>> record = Record(id='NP_055309.2',
            ...                 seq='MRELEAKAT',
            ...                 desc='TNRC6A')
            >>> print(record)
            >NP_055309.2 TNRC6A
            MRELEAKAT
        """
        return self.format(wrap=70).rstrip()

    def __len__(self):
        """Return the length of the sequence.

        Example:
            >>> record = Record(id='NP_055309.2',
            ...                 seq='MRELEAKAT',
            ...                 desc='TNRC6A')
            >>> len(record)
            9
        """
        return len(self.seq)

    def format(self, wrap:int = 70):
        """Returns a formatted FASTA record.

        Args:
            wrap:
                Optional line length to wrap sequence lines (default: 70 
                characters). Use zero (or None) for no wrapping, giving
                a single long line for the sequence.

        Example:
            >>> record = Record(id='NP_055309.2',
            ...                 seq='MRELEAKAT',
            ...                 desc='TNRC6A')
            >>> print(record.format())
            >NP_055309.2 TNRC6A
            MRELEAKAT
            >>> print(record.format(wrap=3))
            >NP_055309.2 TNRC6A
            MRE
            LEA
            KAT
        """
        lst = [self.description, '\n']
        if wrap:
            for i in range(0, len(self.seq), wrap):
                lst.append(f'{self.seq[i : i + wrap]}\n')
        else:
            lst.append(self.seq)
            lst.append('\n')
        return "".join(lst)


def parse(filename: typing.Union[str, pathlib.Path]):
    """Iterates over FASTA records in a file.

    Args:
        filename: A name or path of file containing FASTA sequences.

    Returns:
        A generator of Record objects.
    """
    seqid = None
    desc = None
    seq = []
    with get_open_func(filename) as fh:
        for line in fh:
            if line.startswith('>'):
                if seq:
                    yield Record(seqid, "".join(seq), desc)
                    seq = []
                seqid = line.split()[0][1:]
                desc = line[len(seqid)+1:].strip()
            else:
                seq.append(line.strip())
        if seq:
            yield Record(seqid, "".join(seq), desc)


def to_dict(sequences) -> dict:
    """Turns a generator or list of Record objects into a dictionary.

    This function is not suitable for very large sets of sequences as all the
    SeqRecord objects are held in memory.

    Args:
        sequences: an iterator that returns Record objects, or simply a 
          list of SeqRecord objects.

    Returns:
        A dict mapping sequence id (key) to Record object (value).

    Example:
        >>> import fasta
        >>> record_dict = fasta.to_dict(fasta.parse('test.fasta'))
        >>> print(sorted(record_dict.keys()))
        ['ENO94161.1', 'NP_002433.1', 'sequence']
        >>> print(record_dict['ENO94161.1'].description)
        RRM domain-containing RNA-binding protein
        >>> len(pdict)
        3
    """
    return {record.id: record for record in sequences}

def get_compression_type(
        filename: typing.Union[str, pathlib.Path]
    ) -> typing.Union[str, None]:
    """Returns the compression type of a file based on its first few bytes.

    Args:
        filename: a name or pathlib.Path of a file containing FASTA sequences

    Returns:
        A string representing the compression type of the file, or None
        if the compression type could not be determined.

    Reference:
        http://stackoverflow.com/questions/13044562
    """
    MAGIC_DICT = {b'\x1f\x8b\x08': 'gz',
                  b'\x42\x5a\x68': 'bz2',
                  b'\x50\x4b\x03\x04': 'zip'}
    with open(filename, 'rb') as fh:
        file_start = fh.read(max(len(few_bytes) for few_bytes in MAGIC_DICT))
    for first_bytes, compression_type in MAGIC_DICT.items():
        if file_start.startswith(first_bytes):
            return compression_type
    return None

def get_open_func(filename: typing.Union[str, pathlib.Path]) -> typing.Callable:
    """Returns a function to open a file.

    Raises:
        If a compression type of the file is zstandard and zstandard package is
        not installed, an error is raised.
    """
    compression_type = get_compression_type(filename)
    if not compression_type:
        return open(filename)
    elif compression_type == 'gz':
        return gzip.open(filename, 'rt')
    elif compression_type == 'bz2':
        return bz2.open(filename, 'rt')
    elif compression_type == 'zip':
        with zipfile.ZipFile(filename) as z:
            filename = z.namelist()[0]
            return io.TextIOWrapper(z.open(filename))
    elif compression_type == 'zstandard':
        if HAS_ZSTANDARD:
            return zstandard.open
        else:
            raise ImportError(f'Zstandard is required for: {filename}')