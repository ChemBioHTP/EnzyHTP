from helper import MutaFlag_decode

def test_string():
    x = MutaFlag_decode('XA11Y')
    for item in x:
        try:
            print('for',item,':')
            if item[0].isalpha():
                print(item[0],'is alpha')
            if item[1].isalpha():
                print(item[1],'is alpha')
            if item[2].isdigit():
                print(item[2],'is digit')
            if item[3].isalpha():
                print(item[3],'is alpha')
            print('Pass\n')
        except:
            raise Error('MutaFlag produced by decode function is wrong')

def test_strings():
    x = MutaFlag_decode('AB12C,AB222A')
    for item in x:
        try:
            print('for',item,':')
            if item[0].isalpha():
                print(item[0],'is alpha')
            if item[1].isalpha():
                print(item[1],'is alpha')
            if item[2].isdigit():
                print(item[2],'is digit')
            if item[3].isalpha():
                print(item[3],'is alpha')
            print('Pass\n')
        except:
            raise Error('MutaFlag produced by decode function is wrong')

def test_strings_no_chainid():
    x = MutaFlag_decode('E11F,G22H')
    for item in x:
        try:
            print('for',item,':')
            if item[0].isalpha():
                print(item[0],'is alpha')
            if item[1].isalpha():
                print(item[1],'is alpha')
            if item[2].isdigit():
                print(item[2],'is digit')
            if item[3].isalpha():
                print(item[3],'is alpha')
            print('Pass\n')
        except:
            raise Error('MutaFlag produced by decode function is wrong')

def decode_all_test():
    x = MutaFlag_decode('  XA1 1Y , a b22a,A33333A,  3 dsa dsa 44, 43242343   ')
    print(x)

test_string()

test_strings()

test_strings_no_chainid()

decode_all_test()