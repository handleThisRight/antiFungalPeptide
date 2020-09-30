

## Importing Header Files: 
import string
import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.model_selection import GridSearchCV
from sklearn.svm import SVC
from sklearn.metrics import confusion_matrix, accuracy_score

## UserDefined Functions:
def calcComposition(SequenceSeries, AcidCombination, bond = 0):
    """
	Purpose:
	---
	Function calculates simple composition of Sequences.

	Input Arguments:
	---
	`SequenceSeries`    :   [ pandas.Series ]
		Series of numerous Sequences
    
    `AcidCombination`   :   [ list ]
        list of combinations of find compsition of
    
    `bond`              :   [ integer ]
        definiting the types of pairs (if any) in the combinations
        Default = 0

        0, for AminoAcidComposition (ACC)
        1, for DiPeptideComposition (DPC)
        2, for TriPeptideComposition (TPC)

	Returns:
	---
	`composition`       :   [ dictionary ]
		dictionary of individual compositions of each required combination in the sequence for all sequences.
        Keyed as the index of each sequence.

	"""

    composition = {}
    index = 0
    for seq in SequenceSeries:
        index += 1
        temp = []
        for acid in AcidCombination:
            comp = (seq.count(acid)/(len(seq) - bond)) * 100
            temp.append(comp)
        composition[index] = temp

    return composition




## Main:
if __name__ == "__main__":

    input_file = 'train.csv'
    test_file = 'test.csv'
    output_file = 'outADTN.csv'

    ## Read tinput csv as a pandas dataFrame
    AFPdata = pd.read_csv(input_file)
    
    ## Find the smallest and max length of Sequences: 
    count = 0 
    maxm = minm = 0
    for seq in AFPdata["Sequence"]:
        if count == 0:
            count += 1
            maxm = minm = len(seq)
        else:
            if len(seq) > maxm:
                maxm = len(seq)
            elif len(seq) < minm :
                minm = len(seq)
    
    ## List of all individual AminoAcid:
    AminoAcids = list(string.ascii_uppercase)
    AminoAcids.remove('B')
    AminoAcids.remove('J')
    AminoAcids.remove('O')
    AminoAcids.remove('U')
    AminoAcids.remove('X')
    AminoAcids.remove('Z')

    ## Calculating and adding ACC (Amino Acid Composition) to the dataFrame:
    ACC = calcComposition(AFPdata['Sequence'], AminoAcids, 0)
    X_ACC = pd.DataFrame(ACC).T
    print("\n\n************\nACC Calculated\n************\n\n")

    ## List of all Dipeptide Combinations of AminoAcids:
    dipeptideComb = []
    for a in AminoAcids:
        for b in AminoAcids:
            dipeptideComb.append(a+b)
    
    ## Calculating and adding Dipeptide Composition to the dataFrame:
    DPC = calcComposition(AFPdata['Sequence'], dipeptideComb, 1)
    X_DPC = pd.DataFrame(DPC).T
    print("\n\n************\nDPC Calculated\n************\n\n")
    
    ## List of all Tripeptide Combinations of AminoAcids:
    tripeptideComb = []
    for a in AminoAcids:
        for b in AminoAcids:
            for c in AminoAcids:
                tripeptideComb.append(a+b+c)

    ## Calculating and adding Tripeptide Composition to the dataFrame:
    TPC = calcComposition(AFPdata['Sequence'], tripeptideComb, 2)
    X_TPC = pd.DataFrame(TPC).T
    print("\n\n************\nTPC Calculated\n************\n\n")
    
    ## Dictionary of binary profile of each AminoAcid:
    binaryProfile = {}
    for acid in AminoAcids:
        num = ''
        for i in range(len(AminoAcids)):
            if acid == AminoAcids[i]:
                num += '1'
            else:
                num += '0'
        binaryProfile[acid] = num

    ## Calculating N10C10 bianry profile:
    N10C10 = {}
    index = 0
    for seq in AFPdata['Sequence']:
        index += 1
        temp = []
        if (len(seq) >= 10):
            for i in range(10):
                profile = binaryProfile[seq[i]]
                for j in profile:
                    temp.append(j)
            for i in range(10):
                profile = binaryProfile[seq[(-1*i)]]
                for j in profile: 
                    temp.append(j)
        else:
            count = 0
            for i in range(len(seq)):
                count += 1
                profile = binaryProfile[seq[i]]
                for j in profile:
                    temp.append(j)
            for i in range((10-count)):
                profile = '00000000000000000000'
                for j in profile:
                    temp.append(j)
            count = 0
            for i in range(len(seq)):
                count += 1
                profile = binaryProfile[seq[i]]
                for j in profile:
                    temp.append(j)
            for i in range((10-count)):
                profile = '00000000000000000000'
                for j in profile:
                    temp.append(j)
        N10C10[index] = temp
    X_N10C10 = pd.DataFrame(N10C10).T
    print("\n\n************\nN10C10 Calculated\n************\n\n")

    ######################      MODEL      #################################
    ## Distribute data into input features and output lables:
    X = pd.concat([X_ACC, X_DPC, X_TPC, X_N10C10], axis=1 )
    Y = AFPdata['Lable']
    
    ## Spliting data for external validation:
    X_train, X_test, Y_train, Y_test = train_test_split(X, Y, test_size = 0.25)

    ## creating a parameters dictionary for multiple model settings:
    $parameters = {'C': [5, 7, 10 ], 'gamma': [0.1, 0.01, 0.001]}

    ## Running the model(s):
    #classifier = GridSearchCV(SVC(kernel='rbf', cache_size=4000), parameters, verbose = 1)
    classifier = SVC(C = 10, kernel='rbf', gamma= 0.001)

    ## Fitting the model:
    classifier.fit(X, Y)

    print("\n\n************\nModelFit\n************\n\n")

    ###########################################################################

    ## Reading the test File:
    testData = pd.read_csv(test_file)

    ## Making features:
    test_ACC = pd.DataFrame(calcComposition(testData['Sequence'], AminoAcids, 0)).T
    print("\n\n************\nACC Calculated\n************\n\n")
    test_DPC = pd.DataFrame(calcComposition(testData['Sequence'], dipeptideComb, 1)).T
    print("\n\n************\nDPC Calculated\n************\n\n")
    test_TPC = pd.DataFrame(calcComposition(testData['Sequence'], tripeptideComb, 3)).T
    print("\n\n************\nTPC Calculated\n************\n\n")
    N10C10 = {}
    index = 0
    for seq in testData['Sequence']:
        index += 1
        temp = []
        if (len(seq) >= 10):
            for i in range(10):
                profile = binaryProfile[seq[i]]
                for j in profile:
                    temp.append(j)
            for i in range(10):
                profile = binaryProfile[seq[(-1*i)]]
                for j in profile: 
                    temp.append(j)
        else:
            count = 0
            for i in range(len(seq)):
                count += 1
                profile = binaryProfile[seq[i]]
                for j in profile:
                    temp.append(j)
            for i in range((10-count)):
                profile = '00000000000000000000'
                for j in profile:
                    temp.append(j)
            count = 0
            for i in range(len(seq)):
                count += 1
                profile = binaryProfile[seq[i]]
                for j in profile:
                    temp.append(j)
            for i in range((10-count)):
                profile = '00000000000000000000'
                for j in profile:
                    temp.append(j)
        N10C10[index] = temp
    test_N10C10 = pd.DataFrame(N10C10).T
    print("\n\n************\nN10C10 Calculated\n************\n\n")

    X_test = pd.concat([test_ACC, test_DPC, test_TPC, test_N10C10], axis = 1)

    Y_pred = classifier.predict(X_test)
    Y_pred = pd.Series(Y_pred, name = 'Label')
    print("\n\n************\nPrediction Calculated\n************\n\n")
    outputDataFrame = pd.DataFrame()

    outputDataFrame = pd.concat([testData['ID'], Y_pred], axis = 1)

    outputDataFrame.to_csv(output_file, index=False)

    pass
