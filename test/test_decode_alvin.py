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
    x = MutaFlag_decode('XA11Y,AB22A')
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
    x = MutaFlag_decode('X11Y,A22A')
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

def decode_print():
    x = MutaFlag_decode('XA11Y,AB22A,A33333A,3dsadsa44,43242343')
    print(x)

test_string()

test_strings()

test_strings_no_chainid()

decode_print()