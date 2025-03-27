import os
import numpy as np
from Bio import SeqIO
import shutil
import requests
import pandas as pd
from flask import Flask, render_template, request, send_file, flash, redirect, url_for
from werkzeug.utils import secure_filename
from datetime import datetime

app = Flask(__name__)
app.secret_key = 'supersecretkey'

# Configure upload folder
UPLOAD_FOLDER = 'uploads'
ALLOWED_EXTENSIONS = {'ab1'}
ALLOWED_FOLDER_EXTENSIONS = {'ab1', ''}  # Allow directories

app.config['UPLOAD_FOLDER'] = UPLOAD_FOLDER
os.makedirs(UPLOAD_FOLDER, exist_ok=True)


# Helper function to find AB1 files recursively
def find_ab1_files(path):
    ab1_files = []
    for root, dirs, files in os.walk(path):
        for file in files:
            if file.lower().endswith('.ab1'):
                ab1_files.append(os.path.join(root, file))
    return ab1_files


# Reimplementing functions from the original script
def movingaverage(data_in, window_size):
    window = np.ones(int(window_size)) / float(window_size)
    return np.convolve(data_in, window, 'same')

def process_ab1_file(ab1_file, window_size=10, qual_cutoff=30):
    # Step 1: Read AB1 file and create fasta and qual files
    sample_seq = SeqIO.read(ab1_file, "abi")
    sample_seq.id = sample_seq.name
    
    # Write fasta and qual files
    fasta_file = os.path.join(UPLOAD_FOLDER, sample_seq.id + ".fasta")
    qual_file = os.path.join(UPLOAD_FOLDER, sample_seq.id + ".qual")
    
    SeqIO.write(sample_seq, fasta_file, "fasta")
    SeqIO.write(sample_seq, qual_file, "qual")
    
    # Step 2: Trim sequences based on quality
    sample_qual = SeqIO.read(qual_file, "qual")
    sample_qual_score = sample_qual.letter_annotations["phred_quality"]
    sample_qual_MA = np.array(movingaverage(sample_qual_score, window_size))
    
    if np.max(sample_qual_MA) > qual_cutoff:
        qual_above = list(np.where(sample_qual_MA > qual_cutoff))[0]
        sample_qual_min = np.min(qual_above)
        sample_qual_max = np.max(qual_above)
        
        sample_qual_trim = sample_qual[sample_qual_min:sample_qual_max]
        sample_seq_trim = sample_seq[sample_qual_min:sample_qual_max]
        
        trim_qual_file = os.path.join(UPLOAD_FOLDER, sample_qual.id + ".trim.qual")
        trim_fasta_file = os.path.join(UPLOAD_FOLDER, sample_seq.id + ".trim.fasta")
        
        SeqIO.write(sample_qual_trim, trim_qual_file, "qual")
        SeqIO.write(sample_seq_trim, trim_fasta_file, "fasta")
        
        # Step 3: Create .fa file with the same content as the trimmed .fasta file
        fa_file = os.path.join(UPLOAD_FOLDER, sample_seq.id + ".fa")
        shutil.copy(trim_fasta_file, fa_file)
        
        return {
            "fasta": fasta_file,
            "qual": qual_file,
            "trim_fasta": trim_fasta_file,
            "trim_qual": trim_qual_file,
            "fa": fa_file
        }
    else:
        print(f"The maximum quality score for {sample_qual.id} is below the cutoff ({qual_cutoff})")
        return None
    

# Modified processing function for single file
def process_single_file(filepath, upload_folder, annotation_scheme):
    try:
        base_name = os.path.splitext(os.path.basename(filepath))[0]
        fasta_result = process_ab1_file(filepath)
        
        if not fasta_result:
            return None, "Failed to process AB1 file"
            
        sequences = read_fasta(fasta_result['fa'])
        sequences_frames = gen_frames(sequences)
        proteins = find_prots(sequences_frames)
        best_protein = list(proteins.values())[0]
        annotation_obj = annotate(best_protein, annotation_scheme)
        annotation_result = annotation_obj.retrieve()
        
        # Generate output filenames
        output_files = {
            'output_fasta': f"{base_name}_output_fasta.fa",
            'best_frame': f"{base_name}_best_frame.fa",
            'annotation_result': f"{base_name}_annotation_result.txt",
            'vh_vl_fasta': f"{base_name}_vh_vl_full_length.fa"
        }
        
        # Save individual files
        with open(os.path.join(upload_folder, output_files['output_fasta']), 'w') as f:
            for header, seq in sequences.items():
                f.write(f"{header}\n{seq}\n")
        
        with open(os.path.join(upload_folder, output_files['best_frame']), 'w') as f:
            f.write(f">{base_name}_Best_Protein_Frame\n{best_protein}\n")
        
        with open(os.path.join(upload_folder, output_files['annotation_result']), 'w') as f:
            f.write(f"Annotation Scheme: {annotation_scheme}\n\n")
            f.write("Regions:\n")
            f.write(str(annotation_result[0]) + "\n\n")
            f.write("Number mapping:\n")
            f.write(str(annotation_result[1]) + "\n")
        
        # Excel data collection
        chain_type = 'Heavy' if any(k.startswith('H-') for k in annotation_result[0].keys()) else 'Light'
        cdr1 = annotation_result[0].get(f'{chain_type[0]}-CDR1', '')
        cdr2 = annotation_result[0].get(f'{chain_type[0]}-CDR2', '')
        cdr3 = annotation_result[0].get(f'{chain_type[0]}-CDR3', '')
        full_sequence = ''.join(annotation_result[1].values())
        
        excel_data = {
            'File Name': os.path.basename(filepath),
            'Heavy/Light Chain': chain_type,
            'VH/VL Full Sequence': full_sequence,
            'CDR1': cdr1,
            'CDR2': cdr2,
            'CDR3': cdr3,
            'Processing Time': datetime.now().strftime('%Y-%m-%d %H:%M:%S')
        }
        
        # Save VH/VL FASTA
        with open(os.path.join(upload_folder, output_files['vh_vl_fasta']), 'w') as f:
            f.write(f">{base_name}_VHVL\n{full_sequence}\n")
        
        return excel_data, output_files
        
    except Exception as e:
        return None, str(e)

def read_fasta(fastafile):
    """
    Reads a fasta file and returns a dictionary with header as keys and sequence as values.
    """
    sequences = []
    with open(fastafile, "r") as f:
        for line in f:
            sequences.append(line.rstrip("\n"))
    
    seq_id = [line for line in sequences if line.startswith(">")]
    seq_id_index = [sequences.index(header) for header in seq_id]
    
    seq_dic = {}
    for i in range(len(seq_id_index)):
        if i == (len(seq_id_index) - 1):
            seq_dic[seq_id[i]] = sequences[seq_id_index[i]+1:]
        else:
            seq_dic[seq_id[i]] = sequences[seq_id_index[i]+1:seq_id_index[i+1]]
    
    seq_dic_2 = {}
    for header, lines in seq_dic.items():
        seq_dic_2[header] = "".join(lines)
    
    return seq_dic_2

def swap_dna(dnastring):
    table = {
        'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
        'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
        'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
        'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
        'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
        'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
        'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
        'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
        'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
        'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
        'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
        'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
        'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
        'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
        'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_',
        'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W',
    }
    protein = []
    end = len(dnastring) - (len(dnastring) % 3) - 1
    for i in range(0, end, 3):
        codon = dnastring[i:i+3]
        if codon in table:
            protein.append(table[codon])
        else:
            protein.append("N")
    return "".join(protein)

def rev_seq(seq):
    trans = []
    for i in seq:
        if i == 'A':
            trans.append('T')
        elif i == 'C':
            trans.append('G')
        elif i == 'G':
            trans.append('C')
        elif i == 'T':
            trans.append('A')
        else:
            trans.append(i)
    trans = ''.join(trans)
    seq_rev = trans[::-1]
    return seq_rev

def frame_id(seq):
    """
    Generates the six reading frames for the given sequence.
    """
    frames = {'+1':[],'right':[],'left':[],'-1':[],'+2':[],'+3':[]}
    seq_rev = rev_seq(seq)
    for j in range(3):
        temp = seq[j:]
        temp_rev = seq_rev[j:]
        seq_trans = swap_dna(temp)
        seq_rev_trans = swap_dna(temp_rev)
        if j == 0:
            frames['+1'] = seq_trans
            frames['-1'] = seq_rev_trans
        if j == 1:
            frames['+2'] = seq_trans
            frames['+3'] = seq_rev_trans
        if j == 2:
            frames['right'] = seq_trans
            frames['left'] = seq_rev_trans
    return frames

def gen_frames(dictionary):
    all_dict = {}
    for key, value in dictionary.items():
        all_dict[key] = frame_id(dictionary[key])
    return all_dict

def find_prots(dictionary):
    prots_dict = {}
    for key, frames in dictionary.items():
        poss_protein = []
        for f in frames:
            poss_protein += oframe(frames[f])
        best = ""
        max_len = 0
        for prot in poss_protein:
            if len(prot) > max_len:
                best = prot
                max_len = len(prot)
        prots_dict[key] = best
    return prots_dict

def oframe(amino):
    oframes = []
    for i in range(len(amino)):
        if amino[i] == 'M':
            temp = amino[i:]
            stop = temp.find('_')
            if stop != -1:
                oframes.append(temp[:stop+1])
            else:
                oframes.append(temp)
    return oframes

class annotate:
    """
    Class `annotate`.
    
    Initiator __init__ takes:
      :param aaseq: STRING: A single-letter amino acid sequence (complete VH or VL chain).
      :param scheme: STRING: "kabat", "chothia", "contact", or "imgt" (lowercase).
    """
    def __init__(self, aaseq, scheme):
        self.aaseq = aaseq
        self.scheme = scheme
    
    def __repr__(self):
        return "Annotation of VH or VL sequence using Kabat, Chothia, Contact, or IMGT scheme"
    
    def output(self, chain, lst, regionlst):
        """
        Prints the FR and CDR regions and returns a list of 2 dictionaries.
        """
        self.chain = chain
        self.lst = lst
        self.regionlst = regionlst

        self.regiondict, self.numberdict = {}, {}
        for i in range(0, len(self.lst), 2):
            self.numberdict[self.lst[i]] = self.lst[i+1]
        
        if self.scheme == "kabat":
            print("Annotation in Kabat scheme:")
        elif self.scheme == "chothia":
            print("Annotation in Chothia scheme:")
        elif self.scheme == "contact":
            print("Annotation in Contact scheme:")
        else:
            print("Annotation in IMGT scheme:")
        
        if self.chain == "L":
            print("L-FR1:  ", self.regionlst[0])
            print("L-CDR1: ", self.regionlst[1])
            print("L-FR2:  ", self.regionlst[2])
            print("L-CDR2: ", self.regionlst[3])
            print("L-FR3:  ", self.regionlst[4])
            print("L-CDR3: ", self.regionlst[5])
            print("L-FR4:  ", self.regionlst[6])
            
            for region, seq in zip(["L-FR1", "L-CDR1", "L-FR2", "L-CDR2", "L-FR3", "L-CDR3", "L-FR4"], self.regionlst):
                self.regiondict[region] = seq
            
            return [self.regiondict, self.numberdict]
        else:
            print("H-FR1:  ", self.regionlst[0])
            print("H-CDR1: ", self.regionlst[1])
            print("H-FR2:  ", self.regionlst[2])
            print("H-CDR2: ", self.regionlst[3])
            print("H-FR3:  ", self.regionlst[4])
            print("H-CDR3: ", self.regionlst[5])
            print("H-FR4:  ", self.regionlst[6])
            
            for region, seq in zip(["H-FR1", "H-CDR1", "H-FR2", "H-CDR2", "H-FR3", "H-CDR3", "H-FR4"], self.regionlst):
                self.regiondict[region] = seq
            
            return [self.regiondict, self.numberdict]
    
    def analyze(self, chain, lst):
        """
        Determines the CDR and FR regions based on the numbered sequence.
        """
        self.chain = chain
        self.lst = lst
        if self.chain == "L":
            self.L_FR1, self.L_CDR1, self.L_FR2, self.L_CDR2, self.L_FR3, self.L_CDR3, self.L_FR4 = ["" for i in range(7)]
            try:
                if self.scheme in ["kabat", "chothia"]:
                    self.L_FR1 = "".join([self.lst[i+1] for i in range(0, self.lst.index("L24"), 2)])
                    self.L_CDR1 = "".join([self.lst[i+1] for i in range(self.lst.index("L24"), self.lst.index("L35"), 2)])
                    self.L_FR2 = "".join([self.lst[i+1] for i in range(self.lst.index("L35"), self.lst.index("L50"), 2)])
                    self.L_CDR2 = "".join([self.lst[i+1] for i in range(self.lst.index("L50"), self.lst.index("L57"), 2)])
                    self.L_FR3 = "".join([self.lst[i+1] for i in range(self.lst.index("L57"), self.lst.index("L89"), 2)])
                    self.L_CDR3 = "".join([self.lst[i+1] for i in range(self.lst.index("L89"), self.lst.index("L98"), 2)])
                    self.L_FR4 = "".join([self.lst[i+1] for i in range(self.lst.index("L98"), len(self.lst), 2)])
                
                return [self.L_FR1, self.L_CDR1, self.L_FR2, self.L_CDR2, self.L_FR3, self.L_CDR3, self.L_FR4]
            except ValueError:
                print("Unable to retrieve complete V region. Make sure the sequence has complete V region")
            except:
                print("An error occurred in the `analyze()` method")
        else:
            self.H_FR1, self.H_CDR1, self.H_FR2, self.H_CDR2, self.H_FR3, self.H_CDR3, self.H_FR4 = ["" for i in range(7)]
            try:
                if self.scheme == "kabat":
                    self.H_FR1 = "".join([self.lst[i+1] for i in range(0, self.lst.index("H31"), 2)])
                    self.H_CDR1 = "".join([self.lst[i+1] for i in range(self.lst.index("H31"), self.lst.index("H36"), 2)])
                    self.H_FR2 = "".join([self.lst[i+1] for i in range(self.lst.index("H36"), self.lst.index("H50"), 2)])
                    self.H_CDR2 = "".join([self.lst[i+1] for i in range(self.lst.index("H50"), self.lst.index("H66"), 2)])
                    self.H_FR3 = "".join([self.lst[i+1] for i in range(self.lst.index("H66"), self.lst.index("H95"), 2)])
                    self.H_CDR3 = "".join([self.lst[i+1] for i in range(self.lst.index("H95"), self.lst.index("H103"), 2)])
                    self.H_FR4 = "".join([self.lst[i+1] for i in range(self.lst.index("H103"), len(self.lst), 2)])
                elif self.scheme == "chothia":
                    self.H_FR1 = "".join([self.lst[i+1] for i in range(0, self.lst.index("H26"), 2)])
                    self.H_CDR1 = "".join([self.lst[i+1] for i in range(self.lst.index("H26"), self.lst.index("H33"), 2)])
                    self.H_FR2 = "".join([self.lst[i+1] for i in range(self.lst.index("H33"), self.lst.index("H52"), 2)])
                    self.H_CDR2 = "".join([self.lst[i+1] for i in range(self.lst.index("H52"), self.lst.index("H57"), 2)])
                    self.H_FR3 = "".join([self.lst[i+1] for i in range(self.lst.index("H57"), self.lst.index("H95"), 2)])
                    self.H_CDR3 = "".join([self.lst[i+1] for i in range(self.lst.index("H95"), self.lst.index("H103"), 2)])
                    self.H_FR4 = "".join([self.lst[i+1] for i in range(self.lst.index("H103"), len(self.lst), 2)])
                
                return [self.H_FR1, self.H_CDR1, self.H_FR2, self.H_CDR2, self.H_FR3, self.H_CDR3, self.H_FR4]
            except ValueError:
                print("Unable to retrieve complete V region. Make sure the sequence has complete V region")
            except:
                print("An error occurred in the `analyze()` method")
    
    def retrieve(self):
        """
        Retrieve numbered residues from the Abnum website.
        """
        self.url = "http://www.bioinf.org.uk/abs/abnum/abnum.cgi"
        try:
            if self.scheme not in ["kabat", "chothia", "contact", "imgt"]:
                raise Exception
        except ValueError:
            print("Incorrect scheme mode. Must be one of: kabat, chothia, contact, imgt (in lowercase)")
        else:
            if self.scheme == "kabat":
                self.sche = "-k"
            else:
                self.sche = "-c"
        
        try:
            self.d = {"plain": 1, "scheme": self.sche, "aaseq": self.aaseq}
            self.myPage = requests.get(self.url, params=self.d)
            self.text = self.myPage.text
            self.lst = self.text.split()
                
            if len(self.lst) > 1:
                self.chain = self.lst[0][0]
                self.result = self.output(self.chain, self.lst, self.analyze(self.chain, self.lst))
                return self.result
            else:
                print("No annotation retrieved. Did you enter the complete VH or VL sequence?")
        except:
            print("An error occurred in the `retrieve()` method")

def allowed_file(filename):
    return '.' in filename and \
           filename.rsplit('.', 1)[1].lower() in ALLOWED_EXTENSIONS

@app.route('/', methods=['GET', 'POST'])
def upload_file():
    analysis_results = None
    batch_summary = None
    
    if request.method == 'POST':
        processing_mode = request.form.get('processing_mode', 'single')
        annotation_scheme = request.form.get('annotation_scheme', 'chothia')
        all_excel_data = []
        batch_files = []
        vhvl_entries = []  # Initialize collection for FASTA entries
        
        try:
            if processing_mode == 'single':
                # Single file processing
                file = request.files['file']
                if file and allowed_file(file.filename):
                    filename = secure_filename(file.filename)
                    filepath = os.path.join(app.config['UPLOAD_FOLDER'], filename)
                    file.save(filepath)
                    
                    excel_data, output_files = process_single_file(filepath, app.config['UPLOAD_FOLDER'], annotation_scheme)
                    if excel_data:
                        all_excel_data.append(excel_data)
                        batch_files.append(output_files)
            
            elif processing_mode == 'batch':
                uploaded_files = request.files.getlist('files[]')
                for file in uploaded_files:
                    if file and allowed_file(file.filename):
                        filename = secure_filename(file.filename)
                        filepath = os.path.join(app.config['UPLOAD_FOLDER'], filename)
                        file.save(filepath)
                        
                        if os.path.isdir(filepath):
                            ab1_files = find_ab1_files(filepath)
                            for ab1_file in ab1_files:
                                data, outputs = process_single_file(ab1_file, app.config['UPLOAD_FOLDER'], annotation_scheme)
                                if data:
                                    all_excel_data.append(data)
                                    batch_files.append(outputs)
                                    # Collect FASTA entries
                                    base_name = os.path.splitext(os.path.basename(ab1_file))[0]
                                    vhvl_entries.append({
                                        'base_name': base_name,
                                        'sequence': data['VH/VL Full Sequence']
                                    })
                        else:
                            data, outputs = process_single_file(filepath, app.config['UPLOAD_FOLDER'], annotation_scheme)
                            if data:
                                all_excel_data.append(data)
                                batch_files.append(outputs)
                                # Collect FASTA entries
                                base_name = os.path.splitext(os.path.basename(filepath))[0]
                                vhvl_entries.append({
                                    'base_name': base_name,
                                    'sequence': data['VH/VL Full Sequence']
                                })
            
            # Create output files if any data exists
            if all_excel_data:
                timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
                
                # 1. Create Excel summary
                batch_excel_name = f"batch_summary_{timestamp}.xlsx"
                batch_excel_path = os.path.join(app.config['UPLOAD_FOLDER'], batch_excel_name)
                pd.DataFrame(all_excel_data).to_excel(batch_excel_path, index=False)
                
                # 2. Create consolidated FASTA
                batch_fasta_name = None
                if vhvl_entries:
                    batch_fasta_name = f"batch_vhvl_sequences_{timestamp}.fa"
                    batch_fasta_path = os.path.join(app.config['UPLOAD_FOLDER'], batch_fasta_name)
                    with open(batch_fasta_path, 'w') as f:
                        for entry in vhvl_entries:
                            f.write(f">{entry['base_name']}\n{entry['sequence']}\n")
                
                # Prepare batch summary
                batch_summary = {
                    'total_files': len(all_excel_data),
                    'success_count': len(all_excel_data),
                    'error_count': 0,
                    'excel_file': batch_excel_name,
                    'vhvl_fasta': batch_fasta_name,
                    'output_files': batch_files
                }
            
            analysis_results = {
                'batch_summary': batch_summary,
                'individual_results': all_excel_data[:5]
            }
        
        except Exception as e:
            flash(f'Batch processing error: {str(e)}')
    
    return render_template('upload.html', analysis_results=analysis_results)

@app.route('/download/<filename>')
def download_file(filename):
    return send_file(os.path.join(app.config['UPLOAD_FOLDER'], filename), 
                     as_attachment=True)

if __name__ == '__main__':
    app.run(debug=True)
