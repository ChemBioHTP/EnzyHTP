import sys
import numpy as np
atomlist=[]
coorlist=[]
with open(sys.argv[1]) as f:
    for index,line in enumerate(f):
        if index>1:
           atomlist.append(line.split()[0])
           coorlist.append([float(line.split()[1]),float(line.split()[2]),float(line.split()[3])])
coorlist=np.array(coorlist)
f.close()
natm=len(atomlist)
#CtNidx=[]
CdOidx=[]
CsOidx=[]
CsNidx=[]
CdNidx=[]
PdOidx=[]
PsOidx=[]
CdSidx=[]
CsSidx=[]
#CFidx=[]
#CClidx=[]
#CBridx=[]
#CIidx=[]
for i in range(natm-1):
    for j in range(i+1,natm):
        dis=np.linalg.norm(coorlist[i]-coorlist[j])
        #print(i+1,j+1,atomlist[i],atomlist[j],dis)
        if dis<=1.2: #C-=N (nitrile) identifying
           if atomlist[i]=='C' and atomlist[j]=='N':
              print('Nitrile','C',i+1,'N',j+1)
           elif atomlist[i]=='N' and atomlist[j]=='C':
              print('Nitrile','C',j+1,'N',i+1)
        elif dis<=1.3: #C=O C=N identifying
           if atomlist[i]=='C' and atomlist[j]=='O':
              CdOidx.append([i,j])
           elif atomlist[i]=='O' and atomlist[j]=='C':
              CdOidx.append([j,i])
           elif atomlist[i]=='C' and atomlist[j]=='N':
              CdNidx.append([i,j])
           elif atomlist[i]=='N' and atomlist[j]=='C':
              CdNidx.append([j,i])
        elif dis<=1.335: #C=N C-O identifying
        #!!!!!! Outlier: short C-N in beta-lactam like sub 94,97,102,103,112
           if atomlist[i]=='C' and atomlist[j]=='N':
              CdNidx.append([i,j])
           elif atomlist[i]=='N' and atomlist[j]=='C':
              CdNidx.append([j,i])
           elif atomlist[i]=='C' and atomlist[j]=='O':
              CsOidx.append([i,j])
           elif atomlist[i]=='O' and atomlist[j]=='C':
              CsOidx.append([j,i])
        elif dis<=1.57: #C-O C-N P=O C-F identifying
           if atomlist[i]=='C' and atomlist[j]=='O':
              CsOidx.append([i,j])
           elif atomlist[i]=='O' and atomlist[j]=='C':
              CsOidx.append([j,i])
           elif atomlist[i]=='C' and atomlist[j]=='N':
              CsNidx.append([i,j])
           elif atomlist[i]=='N' and atomlist[j]=='C':
              CsNidx.append([j,i])
           elif atomlist[i]=='P' and atomlist[j]=='O':
              PdOidx.append([i,j])
           elif atomlist[i]=='O' and atomlist[j]=='P':
              PdOidx.append([j,i])
           elif atomlist[i]=='C' and atomlist[j]=='F':
              print('C-F','C',i+1,'F',j+1)
           elif atomlist[i]=='F' and atomlist[j]=='C':
              print('C-F','C',j+1,'F',i+1)
        elif dis<=1.67: #P-O identifying
           if atomlist[i]=='P' and atomlist[j]=='O':
              PsOidx.append([i,j])
           elif atomlist[i]=='O' and atomlist[j]=='P':
              PsOidx.append([j,i])
           elif atomlist[i]=='C' and atomlist[j]=='S':
              CdSidx.append([i,j])
           elif atomlist[i]=='S' and atomlist[j]=='C':
              CdSidx.append([j,i])
        elif dis<=2.25: #C-S C-X (except C-F) identifying
           if atomlist[i]=='C' and atomlist[j]=='S':
              CsSidx.append([i,j])
           elif atomlist[i]=='S' and atomlist[j]=='C':
              CsSidx.append([j,i])
           elif atomlist[i]=='C' and atomlist[j]=='Cl':
              print('C-Cl','C',i+1,'Cl',j+1)
           elif atomlist[i]=='Cl' and atomlist[j]=='C':
              print('C-Cl','C',j+1,'Cl',i+1)
           elif atomlist[i]=='C' and atomlist[j]=='Br':
              print('C-Br','C',i+1,'Br',j+1)
           elif atomlist[i]=='Br' and atomlist[j]=='C':
              print('C-Br','C',j+1,'Br',i+1)
           elif atomlist[i]=='C' and atomlist[j]=='I':
              print('C-I','C',i+1,'I',j+1)
           elif atomlist[i]=='I' and atomlist[j]=='C':
              print('C-I','C',j+1,'I',i+1)
nCdO=len(CdOidx)
#nCN=len(CsNidx)
CdOdel=[]
for i in range(nCdO-1):
    for j in range(i+1,nCdO):
        if CdOidx[i][0]==CdOidx[j][0]: #COO- identifying
           #print('COO-','C',CdOidx[i][0]+1,'O',CdOidx[i][1]+1,'O',CdOidx[j][1]+1)
           CdOdel.append(CdOidx[i])
           CdOdel.append(CdOidx[j])
for i in CdOidx:
    for j in CsOidx:
        if i[0]==j[0]: #Ester/Acid identifying
           if i not in CdOdel:
              CdOdel.append(i)
           for k in CsOidx:
               if j[1]==k[1] and j[0]!=k[0]: #Exclude Acid (COOH)
                  print('Ester(C=O)','C',i[0]+1,'O',i[1]+1)
                  print('Ester(C-O)','C',i[0]+1,'O',j[1]+1)
    for j in CsNidx:
        if i[0]==j[0]: #Amide identifying
           print('Amide(C=O)','C',i[0]+1,'O',i[1]+1)
           print('Amide(C-N)','C',i[0]+1,'N',j[1]+1)
           if i not in CdOdel:
              CdOdel.append(i)
if len(CsSidx)!=0 or len(CdSidx)!=0:
   for i in CdOidx+CdSidx:
       for j in CsOidx+CsSidx:
           if i[0]==j[0]:
              if atomlist[i[1]]=='S' and atomlist[j[1]]=='O': #(C=S)-O
                 print('Ester_S(C=S)','C',i[0]+1,'S',i[1]+1)
                 print('Ester_S(C-O)','C',i[0]+1,'O',j[1]+1)
              elif atomlist[i[1]]=='O' and atomlist[j[1]]=='S': #(C=O)-S
                 print('Ester_S(C=O)','C',i[0]+1,'O',i[1]+1)
                 print('Ester_S(C-S)','C',i[0]+1,'S',j[1]+1)
                 if i not in CdOdel:
                    CdOdel.append(i)
              elif atomlist[i[1]]=='S' and atomlist[j[1]]=='S': #(C=S)-S
                 print('Ester_S(C=S)','C',i[0]+1,'S',i[1]+1)
                 print('Ester_S(C-S)','C',i[0]+1,'S',j[1]+1)
for i in CdOdel:
    CdOidx.remove(i)
for i in CdOidx: #Aldehyde/ketone
    print('Aldehyde/ketone','C',i[0]+1,'O',i[1]+1)
for i in CdNidx:
    for j in CsNidx:
        if i[0]==j[0]: #Amidine identifying
           print('Amidine(C=N)','C',i[0]+1,'N',i[1]+1)
           print('Amidine(C-N)','C',i[0]+1,'N',j[1]+1)
for i in PsOidx: #Phosphate
    for j in CsOidx:
        if i[1]==j[1]:
           for k in PdOidx:
               if i[0]==k[0]:
                  print('Phosphate(P=O)','P',i[0]+1,'O',k[1]+1)
           print('Phosphate(P-O-C)','P',i[0]+1,'O',i[1]+1,'C',j[0]+1)
nPsO=len(PsOidx)
for i in range(nPsO-1):
    for j in range(i+1,nPsO):
        if PsOidx[i][1]==PsOidx[j][1]: #Polyphosphate
           for k in PdOidx:
               if PsOidx[i][0]==k[0] or PsOidx[j][0]==k[0]:
                  print('Phosphate(P=O)','P',k[0]+1,'O',k[1]+1)
           print('Phosphate(P-O-P)','P',PsOidx[i][0]+1,'O',PsOidx[i][1]+1,'P',PsOidx[j][0]+1)
nCsO=len(CsOidx)
for i in range(nCsO-1):
    for j in range(i+1,nCsO):
        if CsOidx[i][1]==CsOidx[j][1]: # Ether
           dis=np.linalg.norm(coorlist[CsOidx[i][0]]-coorlist[CsOidx[j][0]])
           if dis<=1.6: #Epoxide identigying
              print('Epoxide','C',CsOidx[i][0]+1,'O',CsOidx[i][1]+1,'C',CsOidx[j][0]+1)
        elif CsOidx[i][0]==CsOidx[j][0]: # Possibly glycosidic bond
           print('Possible_glycosidic(O-C-O)','O',CsOidx[i][1]+1,'C',CsOidx[i][0]+1,'O',CsOidx[j][1]+1)
if len(CsNidx)!=0 or len(CsSidx)!=0:
   for i in CsOidx:
       for j in CsNidx+CsSidx:
           if i[0]==j[0]:
              if atomlist[j[1]]=='N': # Possibly glycosidic bond
                 print('Possible_glycosidic(O-C-N)','O',i[1]+1,'C',i[0]+1,'N',j[1]+1)
              else: # Possibly glycosidic bond
                 print('Possible_glycosidic(O-C-S)','O',i[1]+1,'C',i[0]+1,'S',j[1]+1)
