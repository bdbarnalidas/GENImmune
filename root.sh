###########################################################################                     Notes
# Author - Barnali Das

# Tested with Bebipred-2.0 present in IEDB analysis resource. But there is an internal server error and we get no response. 
# We need to input a threshold value in order to run most of the tools. The pipeline runs by taking the default threshold values and default window sizes provided by the servers.

# Bebipred - BepiPred predicts the location of linear B-cell epitopes using a combination of a hidden Markov model and a propensity scale method. The residues with scores above the threshold (default value is 0.35) are predicted to be part of an epitope 
#            and marked with "E" in the output table. We have incorporated the Bebipred implementation provided by the IEDB Analysis Resource. Default threshold is 0.350.

# Chou-Fasman - Predicts beta turns to predict antibody epitopes. We have incorporated the Chou-Fasman implementation provided by the IEDB Analysis Resource. Deafult window size is 7. Deafult threshold is the average of all residue scores.

# Emini - Calculation is based on surface accessibility scale. Default window size is 6 and default threshold is 1.0.

# Karplus-Schulz - Deafult threshold is the average of all residue scores. Default window size is 7.

# Kolaskar-Tongaonkar - Deafult threshold is the average of all residue scores. Default window size is 7.

# Parker - Deafult threshold is the average of all residue scores. Default window size is 7.

# Discotope-1.1 - Default threshold is -7.7.

# Discotope-2.0 - Default threshold is -3.7.

# Ellipro - Predicts both linear and discontinuous epitopes. Default minimum score is 0.5. Default maximum distance in Angstroms is 6. 

###########################################################################                     Bash script to execute the entire pipeline  
#!/bin/bash

###########################################################################                     Predicting epitopes from the existing tools directly
function direct_prediction
{
    printf '\n----------------------------------------------------------------------------------- Predicting linear epitopes from protein sequences  ---------------------------------------------------------------------------------------------\n\n';
	printf '\n-------------------------------------------------------------------------------------------------------- Bebipred ------------------------------------------------------------------------------------------------------------------\n';
	printf 'Please type the Swiss-Prot IDs or the protein sequences in fasta format in the file Bebipred/Input.txt. In case of multiple entries, please place the ids and the sequences in separate lines without any blank lines in between them.\n'
    printf 'Waiting for you to provide your inputs\n'
    sleep 5s
    printf 'Have you provided your inputs? Then please do the following:\n'
    sleep 2s
    printf 'Please press 1 if you are entering Swiss-Prot IDs in the input file\n'
    printf 'Please press 2 if you are entering protein sequences in fasta format in the input file\n'
    echo "Your choice is: "
    read choice
    if [ $choice -eq  1 ]; 
    then
        python3.9 -c "import function_library; function_library.run_bebipred_swissprot()"
    else
        python3.9 -c "import function_library; function_library.run_bebipred_sequence()"
    fi

    printf '\n-------------------------------------------------------------------------------------------------------- Chou & Fasman -------------------------------------------------------------------------------------------------------------\n';
    echo "Do you want to run Chou & Fasman for the same input file? If yes, please press Y or else press N.";
    read choice
    if [ "$choice" =  "Y" ]; 
    then 
        cp Bebipred/Input.txt ChouFasman/Input.txt
    else
        printf 'Please type the Swiss-Prot IDs or the protein sequences in fasta format in the file ChouFasman/Input.txt. In case of multiple entries, please place the ids and the sequences in separate lines without any blank lines in between them.\n'
        printf 'Waiting for you to provide your inputs\n'
        sleep 5s
        printf 'Have you provided your inputs? Then please do the following:\n'
        sleep 2s
    fi
    printf 'Please press 1 if the input file contains Swiss-Prot ids\n'
    printf 'Please press 2 if the input file contains protein sequences in fasta format\n'
    echo "Your choice is: "
    read choice
    if [ $choice -eq  1 ]; 
    then
        python3.9 -c "import function_library; function_library.run_choufasman_swissprot()"
    else
        python3.9 -c "import function_library; function_library.run_choufasman_sequence()"
    fi

    printf '\n-------------------------------------------------------------------------------------------------------- Emini ---------------------------------------------------------------------------------------------------------------------\n';
    echo "Do you want to run Emini for the same input file? If yes, please press Y or else press N.";
    read choice
    if [ "$choice" =  "Y" ]; 
    then 
        cp Bebipred/Input.txt Emini/Input.txt
    else
        printf 'Please type the Swiss-Prot IDs or the protein sequences in fasta format in the file Emini/Input.txt. In case of multiple entries, please place the ids and the sequences in separate lines without any blank lines in between them.\n'
        printf 'Waiting for you to provide your inputs\n'
        sleep 5s
        printf 'Have you provided your inputs? Then please do the following:\n'
        sleep 2s
    fi
    printf 'Please press 1 if the input file contains Swiss-Prot ids\n'
    printf 'Please press 2 if the input file contains protein sequences in fasta format\n'
    echo "Your choice is: "
    read choice
    if [ $choice -eq  1 ]; 
    then
        python3.9 -c "import function_library; function_library.run_emini_swissprot()"
    else
        python3.9 -c "import function_library; function_library.run_emini_sequence()"
    fi

    printf '\n-------------------------------------------------------------------------------------------------------- Karplus-Schulz --------------------------------------------------------------------------------------------------------------\n';
    echo "Do you want to run Karplus-Schulz for the same input file? If yes, please press Y or else press N.";
    read choice
    if [ "$choice" =  "Y" ]; 
    then 
        cp Bebipred/Input.txt KarplusSchulz/Input.txt
    else
        printf 'Please type the Swiss-Prot IDs or the protein sequences in fasta format in the file KarplusSchulz/Input.txt. In case of multiple entries, please place the ids and the sequences in separate lines without any blank lines in between them.\n'
        printf 'Waiting for you to provide your inputs\n'
        sleep 5s
        printf 'Have you provided your inputs? Then please do the following:\n'
        sleep 2s
    fi
    printf 'Please press 1 if the input file contains Swiss-Prot ids\n'
    printf 'Please press 2 if the input file contains protein sequences in fasta format\n'
    echo "Your choice is: "
    read choice
    if [ $choice -eq  1 ]; 
    then
        python3.9 -c "import function_library; function_library.run_karplusschulz_swissprot()"
    else
        python3.9 -c "import function_library; function_library.run_karplusschulz_sequence()"
    fi

    printf '\n-------------------------------------------------------------------------------------------------------- Kolaskar-Tongaonkar ---------------------------------------------------------------------------------------------------------\n';
    echo "Do you want to run Kolaskar-Tongaonkar for the same input file? If yes, please press Y or else press N.";
    read choice
    if [ "$choice" =  "Y" ]; 
    then 
        cp Bebipred/Input.txt KolaskarTongaonkar/Input.txt
    else
        printf 'Please type the Swiss-Prot IDs or the protein sequences in fasta format in the file KolaskarTongaonkar/Input.txt. In case of multiple entries, please place the ids and the sequences in separate lines without any blank lines in between them.\n'
        printf 'Waiting for you to provide your inputs\n'
        sleep 5s
        printf 'Have you provided your inputs? Then please do the following:\n'
        sleep 2s
    fi
    printf 'Please press 1 if the input file contains Swiss-Prot ids\n'
    printf 'Please press 2 if the input file contains protein sequences in fasta format\n'
    echo "Your choice is: "
    read choice
    if [ $choice -eq  1 ]; 
    then
        python3.9 -c "import function_library; function_library.run_kolaskartongaonkar_swissprot()"
    else
        python3.9 -c "import function_library; function_library.run_kolaskartongaonkar_sequence()"
    fi

    printf '\n-------------------------------------------------------------------------------------------------------- Parker-----------------------------------------------------------------------------------------------------------------------\n';
    echo "Do you want to run Parker for the same input file? If yes, please press Y or else press N.";
    read choice
    if [ "$choice" =  "Y" ]; 
    then 
        cp Bebipred/Input.txt Parker/Input.txt
    else
        printf 'Please type the Swiss-Prot IDs or the protein sequences in fasta format in the file Parker/Input.txt. In case of multiple entries, please place the ids and the sequences in separate lines without any blank lines in between them.\n'
        printf 'Waiting for you to provide your inputs\n'
        sleep 5s
        printf 'Have you provided your inputs? Then please do the following:\n'
        sleep 2s
    fi
    printf 'Please press 1 if the input file contains Swiss-Prot ids\n'
    printf 'Please press 2 if the input file contains protein sequences in fasta format\n'
    echo "Your choice is: "
    read choice
    if [ $choice -eq  1 ]; 
    then
        python3.9 -c "import function_library; function_library.run_parker_swissprot()"
    else
        python3.9 -c "import function_library; function_library.run_parker_sequence()"
    fi

    printf '\n----------------------------------------------------------------------------------- Predicting discontinuous epitopes from protein structures ----------------------------------------------------------------------------------------\n\n';
    printf '\n-------------------------------------------------------------------------------------------------------- Discotope-1.1 -----------------------------------------------------------------------------------------------------------------\n';
    printf 'Please type in file Discotope-1.1/Input.txt, firstly the PDB ID and then their corresponding chain id separated by a single space character. In case of multiple entries, please place the inputs in separate lines without any blank lines in between them.\n'
    printf 'Waiting for you to provide your inputs\n'
    sleep 5s
    python3.9 -c "import function_library; function_library.run_discotope1()"
    printf '\n-------------------------------------------------------------------------------------------------------- Discotope-2.0 -----------------------------------------------------------------------------------------------------------------\n';
    printf 'Please type in file Discotope-2.0/Input.txt, firstly the PDB ID and then their corresponding chain id separated by a single space character. In case of multiple entries, please place the inputs in separate lines without any blank lines in between them.\n'
    printf 'Waiting for you to provide your inputs\n'
    sleep 5s
    python3.9 -c "import function_library; function_library.run_discotope2()"
    printf '\n---------------------------------------------------------------------------------------------------------- Ellipro -------------------------------------------------------------------------------------------------------------------\n\n';
    printf 'Please type the PDB ids in file Ellipro/Input.txt. In case of multiple entries, please place the inputs in separate lines without any blank lines in between them.\n'
    printf 'Waiting for you to provide your inputs\n'
    sleep 5s
    python3.9 -c "import function_library; function_library.run_ellipro()"
}                         

###########################################################################                     Predicting epitopes from the conserved sequences of MSA
function prediction_from_msa
{
	echo "Inside function 2";
}

###########################################################################                     Main menu
printf 'Please provide your choice from the following:\n'
printf '1. Proceed to epitope prediction directly (Press 1 and enter)\n'
printf '2. MSA and predict epitopes from the conserved sequences (Press 2 and enter)\n'
echo "Your choice is: "
read choice

if [ $choice -eq  1 ]; 
then
    echo "Proceeding to epitope prediction directly";
    direct_prediction;

else
    echo "Proceeding to epitope prediction from the conserved sequences of MSA"
    prediction_from_msa;
fi 

###########################################################################                     Clear input files
