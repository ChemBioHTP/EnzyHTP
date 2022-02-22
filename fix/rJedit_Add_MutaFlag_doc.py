def Add_MutaFlag(self,Flag = 'r', if_U = 0, if_self=0):
        """Returns a label of mutations.

        Reads input of single string or list of strings. Generates a random 
        flag if input is "random" or "r". Appends self.MutaFlags with given
        flag.

        Args:
            Flag: A string or a list of strings. 'r' or 'random' chooses a
              random target residue mutation.

                Example input:

                X: original residue letter. X denotes unknown residue.
                A: chain index.
                11: residue index.
                Y: target residue mutation.

            if_U: include mutations to U (selenocysteine) in random generation 
            if_self: for random generation, compares mutation to original 
              residue.

        Returns:
            A list of strings that represent the original residue letter,
            chain index, residue index, and target residue letter, respectively.
        """

        if type(Flag) == str:
            
            # Random flag generation involves random selection from a
            # given structure. No sorting is applied. 

            if Flag == 'r' or Flag == 'random':
                resi_1 = ''
                resi_2 = ''
                Muta_c_id = ''
                Muta_r_id = ''

                self.get_stru() # Initializes the structure.
                
                chain = choice(self.stru.chains)
                resi = choice(chain.residues)
                Muta_c_id = chain.id
                Muta_r_id = str(resi.id)
                
                if resi.name in Resi_map2:
                    resi_1 = Resi_map2[resi.name]
                else:
                    resi_1 = resi.name
                if if_U:
                    m_Resi_list = Resi_list
                else:
                    m_Resi_list = Resi_list[:-1]
                resi_2 = choice(m_Resi_list)
                if not if_self: # Avoids mutations to original residue.
                    while resi_2 == resi_1:
                        resi_2 = choice(m_Resi_list)

                MutaFlag = (resi_1, Muta_c_id, Muta_r_id ,resi_2)
                self.MutaFlags.append(MutaFlag)

            else:
                MutaFlag = self._read_MutaFlag(Flag)
                self.MutaFlags.append(MutaFlag)

        if type(Flag) == list:
            for i in Flag:
                MutaFlag = self._read_MutaFlag(i)
                self.MutaFlags.append(MutaFlag)

        if Config.debug >= 1:
            print('Current MutaFlags:')
            for flag in self.MutaFlags:
                print(self._build_MutaName(flag))

        label=''
        for flag in self.MutaFlags:
            label=label+'_'+self._build_MutaName(flag) 
        return label
