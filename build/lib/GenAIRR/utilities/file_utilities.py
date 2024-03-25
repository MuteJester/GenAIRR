def parse_fasta(file):
    """Parses a FASTA format file.

    Args:
        file (file): A FASTA formatted file, each sequence is expected to have a
            header line starting with ">".

    Yields:
        generator : A generator object that produces each header and sequence
    """

    header, seq = None, []
    for line in file:
        line = line.rstrip("\n")  # remove newline character
        if line.startswith(">"):  # if line is a FASTA header line
            if header:
                yield header, "".join(seq)  # yield header and sequence
            header, seq = line, []  # set header to contents of line and blank sequence
        else:
            seq.append(line)  # if sequence line, append contents to sequence
    if header:
        yield header, "".join(seq)  # yield header and sequence
