import os
import numpy as np
from Bio import SeqIO
import shutil
import requests
import pandas as pd
from flask import Flask, render_template, request, send_file, flash
from werkzeug.utils import secure_filename
from datetime import datetime

app = Flask(__name__)
app.secret_key = 'supersecretkey'

# Configure upload folder
UPLOAD_FOLDER = 'uploads'
ALLOWED_EXTENSIONS = {'ab1'}

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

def movingaverage(data_in, window_size):
    window = np.ones(int(window_size)) / float(window_size)
    return np.convolve(data_in, window, 'same')

def process_ab1_file(ab1_file, window_size=10, qual_cutoff=20):
    """Reads AB1 and produces trimmed fasta/qual/fa or returns None if below cutoff."""
    sample_seq = SeqIO.read(ab1_file, "abi")
    sample_seq.id = sample_seq.name

    fasta_file = os.path.join(UPLOAD_FOLDER, sample_seq.id + ".fasta")
    qual_file  = os.path.join(UPLOAD_FOLDER, sample_seq.id + ".qual")
    SeqIO.write(sample_seq, fasta_file, "fasta")
    SeqIO.write(sample_seq, qual_file, "qual")

    sample_qual = SeqIO.read(qual_file, "qual")
    scores = np.array(movingaverage(sample_qual.letter_annotations["phred_quality"], window_size))
    if np.max(scores) <= qual_cutoff:
        return None

    above = np.where(scores > qual_cutoff)[0]
    lo, hi = np.min(above), np.max(above)
    trimmed_seq  = sample_seq[lo:hi]
    trimmed_qual = sample_qual[lo:hi]

    trim_fasta = os.path.join(UPLOAD_FOLDER, sample_seq.id + ".trim.fasta")
    trim_qual  = os.path.join(UPLOAD_FOLDER, sample_seq.id + ".trim.qual")
    SeqIO.write(trimmed_seq, trim_fasta, "fasta")
    SeqIO.write(trimmed_qual, trim_qual, "qual")

    fa_copy = os.path.join(UPLOAD_FOLDER, sample_seq.id + ".fa")
    shutil.copy(trim_fasta, fa_copy)

    return {
        "fasta":      fasta_file,
        "qual":       qual_file,
        "trim_fasta": trim_fasta,
        "trim_qual":  trim_qual,
        "fa":         fa_copy
    }

# All downstream steps produce full N/A defaults, then overwrite on success:
def process_single_file(filepath, upload_folder, annotation_scheme):
    base = os.path.splitext(os.path.basename(filepath))[0]
    timestamp = datetime.now().strftime('%Y-%m-%d %H:%M:%S')

    # Prepare "N/A" defaults
    excel_data = {
        'File Name':          os.path.basename(filepath),
        'Heavy/Light Chain':  'N/A',
        'VH/VL Full Sequence':'N/A',
        'CDR1':               'N/A',
        'CDR2':               'N/A',
        'CDR3':               'N/A',
        'Processing Time':    timestamp
    }
    output_files = {
        'output_fasta':       'N/A',
        'best_frame':         'N/A',
        'annotation_result':  'N/A',
        'vh_vl_fasta':        'N/A'
    }

    try:
        # --- Step 1: convert & trim AB1 ---
        fres = process_ab1_file(filepath)
        if not fres:
            raise ValueError("Quality below cutoff")

        # --- Step 2: read & translate & annotate ---
        from Bio import SeqIO as _SQ  # just to be explicit
        sequences = read_fasta(fres['fa'])
        frames    = gen_frames(sequences)
        prots     = find_prots(frames)
        best_prot = list(prots.values())[0]

        annotator = annotate(best_prot, annotation_scheme)
        regions, numbers = annotator.retrieve()
        # if retrieve() failed it would have printed & regions remain undefined

        # --- Step 3: dump all results to files ---
        out_fasta = f"{base}_output_fasta.fa"
        bfasta    = f"{base}_best_frame.fa"
        anres     = f"{base}_annotation_result.txt"
        vhvl      = f"{base}_vh_vl_full_length.fa"

        # 3a: original translated sequences
        with open(os.path.join(upload_folder, out_fasta), 'w') as f:
            for h,s in sequences.items():
                f.write(f"{h}\n{s}\n")

        # 3b: best protein frame
        with open(os.path.join(upload_folder, bfasta), 'w') as f:
            f.write(f">{base}_Best_Frame\n{best_prot}\n")

        # 3c: annotation dump
        with open(os.path.join(upload_folder, anres), 'w') as f:
            f.write(f"Scheme: {annotation_scheme}\n\nRegions:\n{regions}\n\nNumbers:\n{numbers}\n")

        # 3d: VH/VL full sequence
        chain_type = 'Heavy' if any(k.startswith('H-') for k in regions.keys()) else 'Light'
        full_seq   = ''.join(numbers.values())
        with open(os.path.join(upload_folder, vhvl), 'w') as f:
            f.write(f">{base}_VHVL\n{full_seq}\n")

        # --- Step 4: prepare Excel row ---
        excel_data.update({
            'Heavy/Light Chain':   chain_type,
            'VH/VL Full Sequence': full_seq,
            'CDR1':                regions.get(f"{chain_type[0]}-CDR1", ''),
            'CDR2':                regions.get(f"{chain_type[0]}-CDR2", ''),
            'CDR3':                regions.get(f"{chain_type[0]}-CDR3", '')
        })

        # --- Step 5: register output filenames ---
        output_files.update({
            'output_fasta':      out_fasta,
            'best_frame':        bfasta,
            'annotation_result': anres,
            'vh_vl_fasta':       vhvl
        })

    except Exception as e:
        # any failure leaves excel_data & output_files with N/A
        print(f"Processing {filepath} failed:", e)

    return excel_data, output_files

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
            print(self.lst)
                
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
    single_mode = True
    single_result = None
    batch_summary = None
    batch_details = []

    if request.method == 'POST':
        mode = request.form.get('processing_mode', 'single')
        scheme = request.form.get('annotation_scheme', 'chothia')
        all_excel = []
        all_outputs = []

        if mode == 'single':
            single_mode = True
            f = request.files.get('file')
            if f and allowed_file(f.filename):
                fname = secure_filename(f.filename)
                path  = os.path.join(app.config['UPLOAD_FOLDER'], fname)
                f.save(path)
                ed, of = process_single_file(path, app.config['UPLOAD_FOLDER'], scheme)
                # --- CREATE CSV SUMMARY FOR SINGLE FILE ---
                ts = datetime.now().strftime("%Y%m%d_%H%M%S")
                csv_name = f"single_summary_{ts}.csv"
                csv_path = os.path.join(app.config['UPLOAD_FOLDER'], csv_name)
                pd.DataFrame([ed]).to_csv(csv_path, index=False)
                single_result = {
                    'excel': ed,
                    'outs': of,
                    'csv': csv_name
                }

        else:
            single_mode = False
            files = request.files.getlist('files[]')
            for f in files:
                if f and allowed_file(f.filename):
                    fname = secure_filename(f.filename)
                    path  = os.path.join(app.config['UPLOAD_FOLDER'], fname)
                    f.save(path)
                    ed, of = process_single_file(path, app.config['UPLOAD_FOLDER'], scheme)
                    all_excel.append(ed)
                    all_outputs.append(of)

            # build summary
            succ = sum(1 for e in all_excel if e['VH/VL Full Sequence'] != 'N/A')
            err  = len(all_excel) - succ
            timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
            excel_name = f"batch_summary_{timestamp}.xlsx"
            excel_path = os.path.join(app.config['UPLOAD_FOLDER'], excel_name)
            pd.DataFrame(all_excel).to_excel(excel_path, index=False)

            # consolidated fasta (only successes)
            seqs = [ {'base':e['File Name'], 'seq':e['VH/VL Full Sequence']} 
                     for e in all_excel if e['VH/VL Full Sequence']!='N/A' ]
            fasta_name = None
            if seqs:
                fasta_name = f"batch_vhvl_{timestamp}.fa"
                with open(os.path.join(app.config['UPLOAD_FOLDER'], fasta_name),'w') as out:
                    for s in seqs:
                        out.write(f">{s['base']}\n{s['seq']}\n")

            batch_summary = {
                'total_files': len(all_excel),
                'success_count': succ,
                'error_count': err,
                'excel_file': excel_name,
                'vhvl_fasta': fasta_name
            }
            batch_details = list(zip(all_excel, all_outputs))

    return render_template('upload.html',
                           single_mode=single_mode,
                           single_result=single_result,
                           batch_summary=batch_summary,
                           batch_details=batch_details)

@app.route('/download/<filename>')
def download_file(filename):
    return send_file(os.path.join(app.config['UPLOAD_FOLDER'], filename),
                     as_attachment=True)

if __name__ == '__main__':
    app.run(debug=True)
