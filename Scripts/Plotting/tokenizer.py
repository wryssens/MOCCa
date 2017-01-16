################################################################################
#
# Custom tokenizer for spwf tables.
#
def tokenizer(fname):
    with open(fname) as f:
        chunk = []
        for line in f:
            if (line.startswith('*')):
                yield chunk
                chunk = []
                continue
            chunk.append(line)   
