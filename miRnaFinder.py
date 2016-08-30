import sys
import os




def build_combination_dict(mirna, global_dict, mer_len=6, min_distance=2, max_distance=8):

    ## For all distances (min to max)
    for distance in range(min_distance - 1, max_distance):

        ## For all positions in the mirna
        for pos in range(0, len(mirna)):

            ## First part
            first_part = mirna[pos : pos+mer_len]

            ## Second part
            second_part = mirna[pos+distance+mer_len : pos+distance+mer_len+mer_len]

            ## Mer combination
            if len(first_part) == mer_len and len(second_part) == mer_len:

                mer_combi = first_part + "-" + str(distance) + "-" + second_part


                ## Add to global dict
                if mer_combi not in global_dict:        ## Check if it's present
                    global_dict[mer_combi] = list()

                global_dict[mer_combi].append(mirna)


    return global_dict




def print_combinations_for_mirnas(mirnas, global_dict, output_file):

    mirna_combi_dict = dict()

    for key in global_dict.keys():

        ## only take unique mer combinations
        if len(global_dict[key]) == 1:
            mirna = global_dict[key][0]
            combi = key

            #print "mirna: " + str(mirna)
            #print "combi: " + str(combi)

            if mirna not in mirna_combi_dict:
                mirna_combi_dict[mirna] = list()

            mirna_combi_dict[mirna].append(combi)


    no_kmer_found_count = 0
    
    with open(output_file, "w") as of:
    
        for mirna in mirnas:
            of.write(mirna)
            of.write("\n")
            if mirna in mirna_combi_dict:
                of.write(str(len(mirna_combi_dict[mirna])) + "\n")
                of.write(str(mirna_combi_dict[mirna]) + "\n")
            else:
                no_kmer_found_count += 1
                of.write("No kmer found\n")


    print "No kmers found for: " + str(no_kmer_found_count)



def run_program(mirnas, output_file):
    global_dict = dict()

    ## Build the combinations
    for mirna in mirnas:
        global_dict = build_combination_dict(mirna, global_dict)


    ## DEBUG
    #for key in global_dict.keys():
    #    print key + "  :: " + str(len(global_dict[key]))

    ## Print to file
    print_combinations_for_mirnas(mirnas, global_dict, output_file)




def test(output_file):
    mirnas = ["ACGTCGAGTGAGCGAGTGGGTTACGAATAAGCTA"]
    mirnas += ["TTCTCGGAGATCTCGCTAGATAGCTAGCTAGCCTAAGACT"]
    mirnas += ["GGTATCGATAAAGCCTCGCATAGACTCGCCTCGAAAGAAGACTCTC"]
    mirnas += ["TTCCGATAGCGCCTCAGAATCGCGCCTTTTAAGAAGCT"]
    
    run_program(mirnas, output_file)



def main(input_file, output_file):

    ## List of mirnas parsed from input file
    mirnas = list()


    ## Bild listof mirnas based on input file
    with open(input_file) as i_f:
        for line in i_f:
            line = line.strip()
            if line.startswith(">"):
                id = line
                seq = i_f.next()
                mirnas.append(seq.strip())


    ## Run the program
    run_program(mirnas, output_file)




#test()
#test_2()



if __name__ == '__main__':
    print "####################"
    print "## MiRNAMerFinder ##"
    print "####################"
    
    
    if len(sys.argv) != 3:
        print "Please run MiRNAMerFinder with:"
        print "MiRNAMerFinder <input_file> <output_file>"
        sys.exit(1)
        
    input_fasta_file = sys.argv[1]
    output_file = sys.argv[2]
    
    if input_fasta_file == "__test__":
        print "Running a test"
        test(output_file)
        sys.exit()
    
    
    if not input_fasta_file or not os.path.exists(input_fasta_file):
        print "Please provide a valid input FASTA file with miRNA sequences.\n"
        sys.exit(1)
        
        
    
    
    
    main(input_fasta_file, output_file)


    print "\n####################"
    print "##      DONE      ##"
    print "####################"


