import pickle

mutants = [
['NA22K', 'EA24K'],
['NA22R', 'SA29K', 'EA24Q'],
['NA22K', 'SA29K', 'EA24V', 'KA162D', 'RA163F'],
['SA29K', 'EA24Q'],
['NA22R', 'SA29R', 'RA163L'],
['SA29K', 'EA24K', 'KA162E', 'RA163L'],
['NA22K', 'SA29K', 'EA24Q', 'KA162M', 'RA163L'],
['NA22K', 'SA29R', 'EA24Q', 'KA162I', 'RA163L'],
['NA22R', 'EA24K', 'KA162I', 'RA163L'],
['NA22K', 'SA29K', 'EA24V', 'KA162D', 'RA163L'],
['NA22K', 'SA29R', 'EA24R', 'KA162M', 'RA163L'],
['SA29K', 'EA24K', 'KA162L', 'RA163M'],
['NA22K', 'SA29R', 'EA24Q', 'KA162E', 'RA163L'],
['NA22K', 'SA29K', 'EA24R', 'KA162I', 'RA163F'],
['NA22R', 'SA29K', 'EA24K', 'KA162L', 'RA163L'],
]

with open("mutant_list.pickle", "wb") as of:
    pickle.dump(mutants, of)