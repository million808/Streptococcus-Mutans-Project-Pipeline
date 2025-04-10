def separate_fasta(input_file, characterized_file, uncharacterized_file):
    with open(input_file, 'r') as infile, \
         open(characterized_file, 'w') as char_file, \
         open(uncharacterized_file, 'w') as unchar_file:
        
        current_seq = []
        is_characterized = False
        
        for line in infile:
            if line.startswith('>'):
                if current_seq:
                    if is_characterized:
                        char_file.write(''.join(current_seq))
                    else:
                        unchar_file.write(''.join(current_seq))
                    current_seq = []
                
                current_seq.append(line)
                is_characterized = 'uncharacterized' not in line.lower()
            else:
                current_seq.append(line)
        
        if current_seq:
            if is_characterized:
                char_file.write(''.join(current_seq))
            else:
                unchar_file.write(''.join(current_seq))

input_file = 'ua5v7_DE_proteins.fasta'
characterized_file = 'characterizedfasta.txt'
uncharacterized_file = 'uncharacterizedfasta.txt'
separate_fasta(input_file, characterized_file, uncharacterized_file)
