import re  # Import the regular expression module for pattern matching

class Protein:
    '''
    The Protein class parses a PDB file to extract information about a protein.
    It provides:
    lazy attributes (calculated only when accessed):
      - The protein name from the TITLE field of the PDB file. (attribute name .name)
      - Amino acid composition of the entire protein. (attribute name .composition_aa)
      - Amino acid composition per chain (as percentages). (attribute name .chain_aa)
      - Total length of the protein in amino acids. (attribute name .length)
    dunder methods:
      - __str__ : String representation of the Protein object showing file path, name, and length when for example, is printed.
      - __len__ : Define the length of the object as the length of the protein in amino acids. This allows using len(protein_instance) to get the protein length.
    Attributes:
      - _file_path (str): Internal atribute containing the PDB file path to be analyzed.

    '''

    def __init__(self, file_path):        
        self._file_path = file_path  # Store the path to the PDB file in the internal attribute _file_path
     
    def __str__(self):
        '''
         String representation of the Protein object showing file path, name, and length when for example, is printed.
        '''
        return f"Protein object for file: {self._file_path}\nName: {self.name}\nLength: {self.length} AAs" 
    def __len__(self):
        '''
        Define the length of the object as the length of the protein in amino acids. This allows using len(protein_instance) to get the protein length.
        '''
        return sum(self.composition_aa.values())# Sum the counts of all amino acids to get total length

    @property  
    def name(self):
        '''
        Define a property to access the protein name. A @poperty is a built-in decorator that create a methods that compute an attribute. 
        It asure the calculations are perofmed only if the poperty is called, making the creation of the object faster. 
        In addition, once it is called the result is stored as a attibute saving computing time if called again
        '''
        """Extracts the protein name from the TITLE field in the PDB file"""
        with open(self._file_path, "r") as file:  # Open PDB file in read mode
            lines = file.readlines()  # Read all lines into a list
        for line in lines:  # Iterate through each line in the file
            string = line.strip()  # Remove whitespace at start/end of line
            if re.match(r"title", string, re.IGNORECASE):  # Look for line starting with "TITLE"
                name = re.search(r"title\s+(.*)", string, re.IGNORECASE)  # Capture everything after "TITLE"
                break  # Stop the for loop once the title line is found
        return name.group(1)  # Return only the captured protein name. Capture groups are those between () in the regex pattern. Goupr(0) is the whole match, group(1) is the first group, etc.

    @property
    def composition_aa(self):
        """Returns amino acid composition across the whole protein (counts of each AA)"""
        with open(self._file_path, "r") as file:  # Open the PDB file
            lines = file.readlines()  # Read all lines
        prot_AAs = {}  # Dictionary to store amino acid counts
        for line in lines:  # Iterate over each line
            string = line.strip()  # Clean whitespace
            if re.match(r"SEQRES", string, re.IGNORECASE):  # Process only SEQRES lines
                AAs = re.findall(r"\b[A-Z]{3}\b", string)  # Find all 3-letter amino acid codes
                for aa in AAs:  # Loop through found amino acids
                    if aa not in prot_AAs:  # If amino acid not in dictionary
                        prot_AAs[aa] = 1  # Add it with count 1
                    else:
                        prot_AAs[aa] += 1  # Otherwise increment its count
        return prot_AAs  # Return dictionary of amino acid composition

    @property
    def chain_aa(self):
        """Returns amino acid composition per chain (as percentages)"""
        with open(self._file_path, "r") as file:  # Open PDB file
            lines = file.readlines()  # Read all lines
        chains_AAs = {}  # Dictionary to hold chain-wise amino acid counts
        for line in lines:  # Loop through each line
            string = line.strip()  # Remove whitespace
            if re.match(r"SEQRES", string, re.IGNORECASE):  # Look for SEQRES lines
                chains = re.findall(r"\b[A-Z]\b", string)  # Extract chain IDs (single uppercase letters)
                AAs = re.findall(r"\b[A-Z]{3}\b", string)  # Extract 3-letter amino acid codes
                if chains[0] not in chains_AAs:  # If chain not in dictionary
                    chains_AAs[chains[0]] = {}  # Initialize a new dictionary for that chain
                for aa in AAs:  # Loop over amino acids
                    if aa not in chains_AAs[chains[0]]:  # If AA not yet recorded for chain....
                        chains_AAs[chains[0]][aa] = 1  # Add with count 1
                    else: # in the case the AA is already in the dictionary....
                        chains_AAs[chains[0]][aa] += 1  # Increment count

        chain_100 = {}  # Dictionary to hold percentages
        for chain_id, aa_dict in chains_AAs.items():  # Loop through each chain
            total_aa = sum(aa_dict.values())  # Total number of amino acids in chain
            chain_100[chain_id] = {}  # Initialize dictionary for percentages
            for AA, count in aa_dict.items():  # Loop through amino acids in chain
                chain_100[chain_id][AA] = round(float(count / total_aa * 100), 2)  # Convert to percentage, rounded to 2 decimals

        return chain_100  # Return dictionary of percentages per chain

    



# Example usage with p53 PDB file
p53 = Protein("/home/alexander-bontempo/Desktop/master/python/Actividad 1/mubio07_programacion_python_act1/1tup.pdb")

p53.name  # Get the protein name from the TITLE field
p53.composition_aa  # Get amino acid counts for the whole protein
p53.chain_aa  # Get amino acid composition by chain (percentages)
