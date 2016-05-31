import re
repl_dict = {
'Y': '(C|T)', 
'R': '(A|G)',
'W':'(A|T)',
'S':'(G|C)',
'K':'(T|G)',
'M':'(C|A)',
'D':'(A|G|T)',
'V':'(A|C|G)',
'H':'(A|C|T)',
'B':'(C|G|T)',
'N':'(A|C|G|T)'}



def forsub(matchobj):
 return repl_dict[matchobj.group(0)]

def deg_replace(str_toreplace):
 return re.sub("[YRWSKMDVHBN]" , forsub, str_toreplace)




