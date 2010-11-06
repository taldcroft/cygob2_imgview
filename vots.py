import asciitable
import re

def parse_header(header_lines):
    """Parse VOTS header fields from 'header_lines', which should be an
    iterable that returns lines of the VOTS table.  Returns an dict
    of header fields where all but 'description' are in turn a numpy
    recarray table."""
    header = {}
    keywords = ('DESCRIPTION::', 'COOSYS::', 'PARAM::', 'FIELD::')
    key = 'none'
    for line in header_lines:
        line = line.strip()
        if line in keywords:
            key = line[:-2].lower()
            continue
        if key not in header:
            header[key] = []
        header[key].append(line)

    for key, lines in header.items():
        # Flatten description key, otherwise parse lines as a table
        if key == 'description':
            header[key] = '\n'.join(lines)
        else:
            for quotechar in ['"', "'"]:
                try:
                    header[key] = asciitable.read(lines, quotechar=quotechar)
                    break
                except asciitable.InconsistentTableError, error:
                    pass
            else:
                raise asciitable.InconsistentTableError(error)

    return header

def parse_table(lines):
    """Parse VOTS table data from 'lines' where lines is an iterable
    returning lines of ASCII table values, delimited by a character
    that ParseTable understands."""
    headerlines = []
    datalines = []
    for line in lines:
        line = line.strip()
        if line.startswith('##'):
            continue
        if line.startswith('#'):
            headerlines.append(line[1:])
        else:
            datalines.append(re.sub(r'\t', ' ', line))

    header = parse_header(headerlines)
    data = asciitable.read(datalines, Reader=asciitable.NoHeader,
                           names=header['field']['name'])

    return header, data

