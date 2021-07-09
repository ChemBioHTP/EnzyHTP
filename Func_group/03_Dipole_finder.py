from collections import OrderedDict
dipout=open('03_Dipole_out.dat','w')
with open('02_Func_list.txt') as lt: 
     for line in lt:
         sub=line.split()[0]
         func=line.split()[1]
         dipout.write('Substrate '+sub+' '+func+'\n')
         filefunc=open(sub+'.txt','r')
         lines=filefunc.readlines()
         filefunc.close()
         funclines=list(OrderedDict.fromkeys(lines))
         filedip=open(sub+'/LMOdip.txt','r')
         diplines=filedip.readlines()
         filedip.close()
         for i in funclines:
             if func==i.split()[0]:
                atm=[]
                atm.append(i.split()[2]+i.split()[1])
                atm.append(i.split()[4]+i.split()[3])
                if len(i.split())==7:
                   atm.append(i.split()[6]+i.split()[5])
                atmnamelist=[]
                dipavelist=[]
                dipcnt=[]
                for dip in diplines:
                    if len(dip.split())>4:
                       if dip.split()[1]=='(' and dip.split()[3]=='-':
                          temp1=dip.split()[2]
                          temp2=dip.split()[4]
                          if (atm[0]==temp1 and atm[1]==temp2) or (atm[1]==temp1 and atm[0]==temp2):
                             if [temp1,temp2] not in atmnamelist:
                                atmnamelist.append([temp1,temp2])
                                dipavelist.append(float(dip.split()[-1]))
                                dipcnt.append(1)
                             else:
                                idx=atmnamelist.index([temp1,temp2])
                                dipavelist[idx]=dipavelist[idx]+float(dip.split()[-1])
                                dipcnt[idx]=dipcnt[idx]+1
                          if len(atm)==3:
                             if (atm[1]==temp1 and atm[2]==temp2) or (atm[2]==temp1 and atm[1]==temp2):
                                if [temp1,temp2] not in atmnamelist:
                                   atmnamelist.append([temp1,temp2])
                                   dipavelist.append(float(dip.split()[-1]))
                                   dipcnt.append(1)
                                else:
                                   idx=atmnamelist.index([temp1,temp2])
                                   dipavelist[idx]=dipavelist[idx]+float(dip.split()[-1])
                                   dipcnt[idx]=dipcnt[idx]+1
                for j in range(len(atmnamelist)):
                    dipout.write(' '+atmnamelist[j][0]+' '+atmnamelist[j][1]+' '+str(dipavelist[j]/dipcnt[j])+'\n')
         dipout.write('\n')
dipout.close()
