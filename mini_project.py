#Python program to analyse pdb files

def information():
    string=''
    aminocountlist = []
    aaseq =[]
    chains=[]
    helixcount =''
    sheetstr = ''
    allaaseq={}
    chain_type = []
    output = []
    translatedallseq={}
    amino_dictionary={"ASN":"N","MET":"M","ASP":"D","ASX":"B","LYS":"K","CYS":"C","GLU":"E","GLN":"Q","GLX":"Z","GLY":"G","HIS":"H","PRO":"P","ILE":"I","LEU":"L","PHE":"F","SER":"S","THR":"T","TRP":"W","ARG":"R","ARG":"R","TYR":"Y","VAL":"V","ALA":"A"}
    
    with open(filename) as file:
        for line in file.readlines():
            strfile=line.split()
            colname = strfile[0]        
            if colname =='TITLE':
                head=strfile[1:]
                if head[0].isalpha():
                    head1=' '.join(head[0:])
                    string+=head1
                if head[0].isdigit():
                    head2=' '.join(head[1:])
                    string+=head2
                
            if colname =='SEQRES':
                aminocounts=f'{strfile[2]}:{strfile[3]}'
                aminocountlist.append(aminocounts)
                chain_aminonum_list={strfile[2]:strfile[4:]}
                aaseq.append(chain_aminonum_list)
                for linelist in strfile[2]:
                     chains.append(linelist)
                        
            if colname == 'HELIX':
                chain = strfile[4]
                helixcount+=chain
                
            if colname == 'SHEET':
                sheet = strfile[5]
                sheetstr += sheet
                
    
    for chain_amino in aaseq:
        for chain,seq in chain_amino.items():
            if chain in allaaseq.keys():
                allaaseq[chain]+=seq
            else:
                allaaseq[chain]=seq

    translatedallseq={k:''.join([amino_dictionary.get(v) for v in v]) for k,v in allaaseq.items()}

    chainset=set(chains)
    chainlist=list(chainset)
    chainlist.sort()

    aminocountlist = set(aminocountlist)
    aminocountlist = list(aminocountlist)
    aminocountlist.sort()

    helixdict={}
    for chain in chains:
        for helix in helixcount:

            if chain== helix:
                helixdict[chain]=helixcount.count(helix)
            elif chain not in helixcount:
                helixdict[chain]=0

    sheetdict={}
    for chain in chains:
        for sheet in sheetstr:
            if chain== sheet:
                sheetdict[chain]=sheetstr.count(sheet)
            elif chain not in sheetstr:
                sheetdict[chain]=0

    aminonodict={}

    for chain in chainlist:
        for num in aminocountlist:
                if chain==num[0]:
                    aminonodict[chain]=num[2:]

    sorted(helixdict.keys())
    print(f'PDB File: {filename}')
    print('TITLE: %s '%string)
    for chain in chainlist:
        print(f" -Chain {chain}")
        print('%3s Number of amino acids: %s' %('',aminonodict.get(chain)))
        print('%3s Number of helix:%8s%s' %('',(helixdict.get(chain)),''))
        print('%3s Number of sheets:%7s%s' %('',(sheetdict[chain]),''))

        if len(translatedallseq.get(chain)) > 50:
            for word in range (0, len(translatedallseq.get(chain)),50):
                print('%3s Sequences:  '%(''),(translatedallseq.get(chain)[0:50]))
                break
            for word in range (0, len(translatedallseq.get(chain)),50):
                print('%14s  '%(''),(translatedallseq.get(chain)[word:word+50]))

        else:
            print('%3s Sequences:  '%(''),(translatedallseq.get(chain)))


            
def histogram():
    quit=''
    
    while quit !='q':
    
        print('Choose an option to order by or q to quit:')
        print('%3s number of amino acids - ascending  (an)' %'')
        print('%3s number of amino acids - descending (dn)' %'')
        print('%3s alphabetically - ascending         (aa)' %'')
        print('%3s alphabetically - descending        (da)' %'')
        sortby=input('order by:')
        aaseqmainstr=[]
        allaaseqlist=[]
        amino={}
        aminolist=[]
        amino_dict={}
        for line in open(filename).readlines():
            listlines= line.split()
            colname=listlines[0]
            if colname == 'SEQRES':
                chain=listlines[4:]
                aaseqmainstr.append(' '.join(chain).split())

        for i in aaseqmainstr:
            allaaseqlist+=i
        amino=set(allaaseqlist)
        aminolist=list(amino)

        for i in aminolist:
            amino_dict[i]=allaaseqlist.count(i)

        if sortby == 'aa' :
            aminolist.sort()
            for i in aminolist:   
                print("{:4s} {:1s} {:2d} {:2s} {:2} {:2}".format(i,'(',(allaaseqlist.count(i)),')',':',('*'*allaaseqlist.count(i))))

        elif sortby == 'da' :
            aminolist.sort(reverse=True)
            for i in aminolist:   
                print("{:4s} {:1s} {:2d} {:2s} {:2} {:2}".format(i,'(',(allaaseqlist.count(i)),')',':',('*'*allaaseqlist.count(i))))

        elif sortby == 'an':
            amino_dict_sorted=sorted(amino_dict.items(), key=lambda x: x[1])
            for j in amino_dict_sorted:
                print("{:4s} {:1s} {:2d} {:2s} {:2s} {:2s}".format((j[0]),'(',(j[1]),')',':',('*'*allaaseqlist.count(j[0]))))

        elif sortby == 'dn' :
            amino_dict_sorted=sorted(amino_dict.items(), key=lambda x: x[1],reverse=True)
            for j in amino_dict_sorted:
                print("{:4s} {:1s} {:2d} {:2s} {:2s} {:2s}".format((j[0]),'(',(j[1]),')',':',('*'*allaaseqlist.count(j[0]))))
        elif sortby in ('AN','DA','AA','DN'):
            print("Case sensitive options !!!!")
        elif sortby == 'q':
            quit+=sortby
            
            break
        
        else:
            print("Invalid option")


            
            
            

def secondary():
    helix_num =''
    sheetstring = ''
    chainlist=[]
    chain_type = []
    output = []
    chainaacountlist = []
    chainaalist =[]
    amino_dict={"ASN":"N","MET":"M","ASP":"D","ASX":"B","LYS":"K","CYS":"C","GLU":"E","GLN":"Q","GLX":"Z","GLY":"G","HIS":"H","PRO":"P","ILE":"I","LEU":"L","PHE":"F","SER":"S","THR":"T","TRP":"W","ARG":"R","ARG":"R","TYR":"Y","VAL":"V","ALA":"A"}
    #Reading file, each line as a list to easily access file items
    try:
        with open(filename) as file:
            for line in file.readlines():
                listlines=line.split()
                colname = listlines[0]
                if colname == 'HEADER':
                    header = listlines[4]
        #Creates a list of concatenated chain to its amino acid counts (chainaaminocounts) and to its amino acid residues(chainaalist)
                if colname =='SEQRES':
                    chainaacounts=f'{listlines[2]}:{listlines[3]}'
                    chainaacountlist.append(chainaacounts)
                    chainamino={listlines[2]:listlines[4:]}
                    chainaalist.append(chainamino)

                    for i in listlines[2]:
                        chainlist.append(i)


                if colname == 'HELIX':
                    chain = listlines[4]
                    helix_num+=chain


                if colname == 'SHEET':
                    sheets = listlines[5]
                    sheetstring += sheets
            chainaadict={}
            for i in  chainaalist:
                for chain,amino in  i.items():
                    if chain in  chainaadict.keys():
                         chainaadict[chain]+=amino
                    else:
                         chainaadict[chain]=amino
            allaaseqtransdict={chain:''.join([amino_dict.get(amino) for amino in amino]) for chain,amino in  chainaadict.items()}
            chain=set(chainlist)
            chain=list(chain)
            chain.sort()
            chainaacountlist = set(chainaacountlist)
            chainaacountlist = list(chainaacountlist)
            chainaacountlist.sort()
            helixdict={}
            for a in chainlist:
                for i in helix_num:
                    if a== i:
                        helixdict[a]=helix_num.count(i)
                    elif a not in helix_num:
                        helixdict[a]=0
            sheetdict={}
            for a in chainlist:
                for i in sheetstring:
                    if a== i:
                        sheetdict[a]=sheetstring.count(i)
                    elif a not in sheetstring:
                        sheetdict[a]=0
            aanumdict={}
            for r in chain:
                for i in chainaacountlist:
                    if r==i[0]:
                        aanumdict[r]=i[2:]
                        
            print(f'Secondary Structure of the PDB id {header}:')
            sorted(helixdict.keys())

            x_option="-"
            open_file = []
            with open(filename) as listlines:
                open_file = listlines.readlines()
            for i in chain:
                print(f"Chain {i}:")
                print(f'(1)')
                symbol_list=list('-'*len(allaaseqtransdict.get(i)))
                symbol_list1=list('-'*len(allaaseqtransdict.get(i)))
                for line in open_file:
                    if line.startswith('HELIX'):
                        if i == line.split()[4]:
                            helix_begin_pos=int(line.split()[5])
                            helix_end_pos=int(line.split()[8])
                            helix_position=line.split()[1]
                            for x, v in enumerate(symbol_list,start=1):
                                if x >= helix_begin_pos and x <= helix_end_pos:
                                    symbol_list[x-1] = '/'                               
                                if x == helix_begin_pos:
                                    if len(helix_position)==1:
                                        symbol_list1[x-1]=helix_position
                                    elif len(helix_position)==2:
                                        symbol_list1[x-1:x+1]=helix_position 
                    if line.startswith('SHEET'):
                        if i == line.split()[5]:
                            begin_sheet_pos=int(line.split()[6])
                            end_sheet_pos=int(line.split()[9])
                            wx = line.split()[1]
                            xy = line.split()[2]
                            sheet_position = wx + xy
                            for a, v in enumerate(symbol_list,start=1):
                                if a >= begin_sheet_pos and a <= end_sheet_pos:
                                    symbol_list[a-1] = '|'
                                if a == begin_sheet_pos:
                                    if sheet_position==1:
                                        symbol_list1[a-1]=sheet_position
                                    else:
                                        symbol_list1[a-1:a+1]=sheet_position
                symbols = ''.join(symbol_list)
                symbolstr = ''.join(symbol_list1)
                if len(allaaseqtransdict.get(i))>80:


                    for k in range(0,len(allaaseqtransdict.get(i)),80):
                        print(allaaseqtransdict.get(i)[k:k+80])
                        print(symbols[k:k+80])
                        print(symbolstr[k:k+80].replace('-', ' '))
                    print(f"({len(allaaseqtransdict.get(i))})" '\n')


                else:
                    print(allaaseqtransdict.get(i))
                    print(symbols)
                    print(symbolstr.replace('-',' '))
                    print(f"({len(allaaseqtransdict.get(i))})" '\n') 
    except:
        print("Corrupted pdb file")
        
            
            

filename='None'
def menu():
    print("*" * 80)
    print("* PDB FILE ANALYZER %60s" %'*')
    print("*" * 80)
    print("* Select an option from below:%50s" %'*')
    print("*%79s" %'*')
    print("* %5s 1) Open a PDB File %22s (O) %26s" %('','','*'))
    print("* %5s 2) Information %26s (I) %26s"  %('','','*'))
    print("* %5s 3) Show histogram of amino acids %8s (H) %26s"  %('','','*'))
    print("* %5s 4) Display Secondary Structure %10s (S) %26s"  %('','','*'))
    print("* %5s 5) Export PDB File %22s (X) %26s"  %('','','*'))
    print("* %5s 6) Exit  %32s (Q) %26s"  %('','','*'))
    print("* %78s" %'*')
    print(f"* {('Current PDB:'+ filename).rjust(76)} * ")
    print("*" * 80)
    return input(":")

option=menu()
memory=''
options=['i','h','s','q','2','3','4','6',]

while option != 'E':
    
    if option.lower() != 'q':
        if option.lower() in ('o','1',):
            path=input("Enter a Valid PATH for a PDB File:")
            import os
            if os.path.exists(path):
                filename=os.path.basename(path)
                
                if filename.endswith('.pdb'):
                    try:
                        with open(filename) as file:
                            line = file.readline().replace('\n', '')
                    except FileNotFoundError:
                        print('No such File or Directory')

                    if line.startswith('HEADER') and len(line) == 80:

                        memory+=filename
                        print(f'The File {filename} has been successfully loaded.')

                        if len(memory)!=0:

                            while option!='E':
                                print("Enter analysis option")
                                choice=menu()
                                if choice.lower() in('i','2'):
                                    information()
                                    
                                elif choice.lower() in ('h','3'):
                                    histogram()
                                    
                                elif choice.lower() in ('s','4'):
                                    secondary()
                                elif choice.lower() in ('x','5'):
                                    print("Option deprecated !!!")
                                    
                                elif choice.lower() in ('q','6'):
                                    option2 =input("Do you want to exit(E) or do you want go back to the menu (M) ")
                                    if option2 == 'M':
                                        pass
                                        
                                    if option2 == 'E':
                                        option='E'
                                        
                                elif choice =='1':

                                    print(f'The file ({filename}) already exists !!.Enter R to replace it or E to exit')
                                    replace=input(":")
                                    if replace.lower()=='r':
                                        filename='None'
                                        break
                                        
                                    elif replace.lower()=='E':
                                        break
                                        
                                        
                                    else:
                                        print("Invalid option")
                                        continue
                                else:
                                    print("Invalid analysis option")                                        
                    else:
                        print("File integrity Failed !!!!")
                else:
                    print("File Format Not Supported. Only .pdb files !!! \nChoose 1 or o to load a file or E to quit") 
                    filename='None'
                    option=menu()
            else:
                print("File path does not exist,Choose 1 or o to load a file or E to quit")
                option=menu()
        elif option in options:
            print("No file loaded, \nChoose 1 or o to load a file or E to quit")
            option=menu()
        else:
            print("Invalid option for loading file \nChoose 1 or o to load a file or E to quit")
            option=menu()
    else:
        option =input("Do you want to exit(E) or do you want go back to the main menu (M)")

        if option == 'M':
            pass
            
        else:
            break



