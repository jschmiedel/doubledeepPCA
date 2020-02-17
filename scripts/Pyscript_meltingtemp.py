from Bio.SeqUtils import MeltingTemp as mt
# from Bio.Seq import Seq
import datatable as dt
PDZseq = 'GGAGGTGGAGCTAGCCCGAGGCGAATTGTGATCCACCGGGGCTCCACGGGCCTGGGCTTCAACATCGTGGGTGGCGAGGACGGTGAAGGCATCTTCATCTCCTTTATCCTGGCCGGGGGCCCTGCAGACCTCAGTGGGGAGCTGCGGAAGGGGGACCAGATCCTGTCGGTCAACGGTGTGGACCTCCGAAATGCCAGCCATGAGCAGGCTGCCATTGCCCTGAAGAATGCGGGTCAGACGGTCACGATCATCGCTCAGTATAAACCATAAAAGCTTCGCAGG'


DT = dt.cbind(dt.Frame(pos=range(1, int((len(PDZseq) - 30) / 3))),
    dt.Frame(seq=['A']),
    dt.Frame(Tm_Wallace=[float(0)]),
    dt.Frame(Tm_GC=[float(0)]),
    dt.Frame(Tm_NN=[float(0)]))


for idx in range(1, int((len(PDZseq) - 30) / 3)):
    idx_seq = PDZseq[((idx - 1) * 3 + 1):(30 + idx * 3 + 1)]
    DT[idx - 1, 'seq'] = idx_seq
    DT[idx - 1, 'Tm_Wallace'] = mt.Tm_Wallace(idx_seq)
    DT[idx - 1, 'Tm_GC'] = mt.Tm_GC(idx_seq)
    DT[idx - 1, 'Tm_NN'] = mt.Tm_NN(idx_seq)

DT.to_csv("dataset/PDZ/PDZ_NM_meltingTemps.csv")
