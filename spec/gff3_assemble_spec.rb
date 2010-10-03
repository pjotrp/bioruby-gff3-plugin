# RSpec for BioRuby-GFF3-Plugin. Run with something like:
#
#   ruby -I ../bioruby/lib/ ~/.gems/bin/spec spec/gff3_assemble_spec.rb 
#
# Copyright (C) 2010 Pjotr Prins <pjotr.prins@thebird.nl>
#
$: << "../lib"

require 'bio/db/gff/gffdb'

include Bio::GFFbrowser

FASTAFILE="test/data/gff/MhA1_Contig1133.fa"  
GFF3FILE="test/data/gff/MhA1_Contig1133.gff3"

describe GFFdb, "Assemble CDS" do
  before :all do 
    gffdb = Bio::GFFbrowser::GFFdb.new(GFF3FILE, :fasta_filename => FASTAFILE)
    @gff = gffdb.assembler
    @gff.parse
    @contigsequence = @gff.sequencelist["MhA1_Contig1133"]
    @componentlist = {}
    @cdslist = {}
    @gff.each_CDS do | id, reclist, component |
      @componentlist[id] = component
      @cdslist[id] = reclist
    end
  end

  it "should have the single contig" do 
    @gff.sequencelist.size.should == 1
    @gff.sequencelist["MhA1_Contig1133"].should_not == nil
    @gff.sequencelist["MhA1_Contig1133"].size.should == 33905
  end
  it "should have a container component" do
    component = @componentlist['cds:MhA1_Contig1133.frz3.gene4']
    component.start.should == 7838
    component.end.should == 8740
  end
  it "should have CDS 7838:7980" do
    recs = @cdslist['cds:MhA1_Contig1133.frz3.gene4']
    cds0 = recs[0]
    cds0.start.should == 7838
    cds0.end.should == 7980
    cds0.frame.should == 0
    cds0.seqname.should == 'MhA1_Contig1133'
  end
  it "should have CDS 8065:8308" do
    recs = @cdslist['cds:MhA1_Contig1133.frz3.gene4']
    cds1 = recs[1]
    cds1.start.should == 8065
    cds1.end.should == 8308
    cds1.frame.should == 1
    cds1.strand.should == '+'
    cds1.seqname.should == 'MhA1_Contig1133'
  end
  # From Wormbase website http://www.wormbase.org/db/gb2/gbrowse/m_hapla/?name=MhA1_Contig1133%3A7838..8740
  # >MhA1_Contig1133:7838..8740
  # atgcgtcctttaacagatgaagaaactgaaaagtttttcaaaaaactttcaaattatatt
  # ggtgacaatattaaacttttattggaaagagaagatggagaatatgtttttcgtttacat
  # aaagacagagtttattattgcaggtttttttaaaattattttatatttaaattaggtctc
  # aatctttataggggattttgtttttgttatttttttttggtttttag>tgaaaaattaatg
  # cgacaagcagcatgtattggacgtaaacaattgggatcttttggaacttgtttgggtaaa
  # ttcacaaaaggagggtctttctttcttcatataacatcattggattatttggcaccttat
  # gctttagcaaaaatttggttaaaaccacaagctgaacaacaatttttatatggaaataat
  # attgttaaatctggtgttggaagaatgagtgaagggattgaagaaaaacaagtaaatatt
  # taattattttttttaaaatggattcctttacttctcaattaaatattaaaagcatatctg
  # tagaagaggttatttatctttaaatcgaaatatacaggaataaataaaaatttaagaaat
  # cataatttagaattctttttctggttatgttagattatttttaaatttttttgtaatttt
  # tttttcgtaatttttttatgagcaaatcccttctctcttaaatattttaataaaaatcta
  # attttataaattataattattttttagggtattattatttataatatgtcagatttacca
  # ttgggttttggagtggctgcaaagggaacattatcttgtagaaaagtagatcctacagct
  # ttagttgttttacatcaatcagatttgggtgaatatattcgaaatgaagagggattaatt

  it "should translate CDS 7838:7980 (in frame 0, + strand)" do
    recs = @cdslist['cds:MhA1_Contig1133.frz3.gene4']
    component = @componentlist['cds:MhA1_Contig1133.frz3.gene4']
    cds0 = recs[0]
    cds0.seqname.should == 'MhA1_Contig1133'
    seq = @gff.assemble(@contigsequence,component.start,[cds0])
    seq.should == "ATGCGTCCTTTAACAGATGAAGAAACTGAAAAGTTTTTCAAAAAACTTTCAAATTATATTGGTGACAATATTAAACTTTTATTGGAAAGAGAAGATGGAGAATATGTTTTTCGTTTACATAAAGACAGAGTTTATTATTGCAGG"
    aaseq = @gff.assembleAA(@contigsequence,component.start,[cds0])
    aaseq.should == "MRPLTDEETEKFFKKLSNYIGDNIKLLLEREDGEYVFRLHKDRVYYCR"
  end
  it "should translate CDS 8065:8308 (in frame 1, + strand)" do
    recs = @cdslist['cds:MhA1_Contig1133.frz3.gene4']
    component = @componentlist['cds:MhA1_Contig1133.frz3.gene4']
    cds1 = recs[1]
    seq = @gff.assemble(@contigsequence,component.start,[cds1], :codonize => false)
    seq.should == "TGAAAAATTAATGCGACAAGCAGCATGTATTGGACGTAAACAATTGGGATCTTTTGGAACTTGTTTGGGTAAATTCACAAAAGGAGGGTCTTTCTTTCTTCATATAACATCATTGGATTATTTGGCACCTTATGCTTTAGCAAAAATTTGGTTAAAACCACAAGCTGAACAACAATTTTTATATGGAAATAATATTGTTAAATCTGGTGTTGGAAGAATGAGTGAAGGGATTGAAGAAAAACAAG"
    seq = @gff.assemble(@contigsequence,component.start,[cds1], :codonize => true)
    seq.should == "GAAAAATTAATGCGACAAGCAGCATGTATTGGACGTAAACAATTGGGATCTTTTGGAACTTGTTTGGGTAAATTCACAAAAGGAGGGTCTTTCTTTCTTCATATAACATCATTGGATTATTTGGCACCTTATGCTTTAGCAAAAATTTGGTTAAAACCACAAGCTGAACAACAATTTTTATATGGAAATAATATTGTTAAATCTGGTGTTGGAAGAATGAGTGAAGGGATTGAAGAAAAACAA"
    aaseq = @gff.assembleAA(@contigsequence,component.start,[cds1])
    # note it should handle the frame shift and direction!
    aaseq.should == "EKLMRQAACIGRKQLGSFGTCLGKFTKGGSFFLHITSLDYLAPYALAKIWLKPQAEQQFLYGNNIVKSGVGRMSEGIEEKQ"
    # we are going to force a change of direction
    # cds1rev = cds1
    # cds1rev.strand = '-'
    # seq = @gff.assemble(@contigsequence,component.start,[cds1rev])
    # seq.should == "ACAAAAAGAAGTTAGGGAAGTGAGTAAGAAGGTTGTGGTCTAAATTGTTATAATAAAGGTATATTTTTAACAACAAGTCGAACACCAAAATTGGTTTAAAAACGATTTCGTATTCCACGGTTTATTAGGTTACTACAATATACTTCTTTCTTTCTGGGAGGAAAACACTTAAATGGGTTTGTTCAAGGTTTTCTAGGGTTAACAAATGCAGGTTATGTACGACGAACAGCGTAATTAAAAAGT"
    # aaseq = @gff.assembleAA(@contigsequence,component.start,[cds1rev])
    # aaseq.should == "EKLMRQAACIGRKQLGSFGTCLGKFTKGGSFFLHITSLDYLAPYALAKIWLKPQAEQQFLYGNNIVKSGVGRMSEGIEEKQ"
  end
  it "should assemble 3 CDSs for MhA1_Contig1133.frz3.gene4"
  # > class=Sequence position=MhA1_Contig1133:27463..29904 (- strand)
  # ATGGACCATC ATGCATTGGT GGAGGAATTA CCAGAAATTG AAAAATTAAC TCCTCAAGAA CGTATTGCAT TAGCTAGAGA
  # ACGCCGTGCT GAACAACTTC GACAGAATGC TGCACGGGAG GCTCAATTGC CAATGCCTGC ACAGCGCCGG CCTCGTCTTC
  # GATTTACACC AGATGTTGCT TTACTTGAGG CAACAGTTAG GGGTGATACC CAAGAAGGTT ATACATAAAG ATTATTGATT
  # TTAAATGAAT TTATTTATTT TTTAGTTGAA AGACTTTTAA TGGAAGGTGT CAATGCTGAT TCACATAATG AGGATGGATT
  # AACACCTTTA CATCAGGCAA AAACCAAATT AATTTTTTTA AATTTATTTT TAGTGTGCCA TTGACAATAA TGAAAGAATT
  # GTTCGTCTTC TGCTTAGGTA CGGAGCTTGT GTTAATGCCA AAGACACTGA ACTTTGGACA CCATTGCACG CAGCTGCATG
  # TTGTGCTTAT ATTGATATTG TTCGATTGCT TATTGCACAG TTAGTTTTTT TTTAATTTTT TTTTTAAATA AATTTCTTAA
  # GTTTTACAGA AATATTTATT TTAAACAAAC GGGACTTCCT TTTAAATTTT TTGTATTTTT AATCTTTACG TATTTTCATT
  # TAATAATTAA TTCGTCTTCT AAAAGTTCGT AAGTTTTGTG GTTTAGTTTA ATGGGTAAAC ATCCAGTTTT TAGGTCATCG
  # ATTTTTATTT TTGCGTCATA TTTTATCGAA AACTTCTTTC ATATTAAAAA TTTCTTTTTA AGCAACGCAG ATTTACTAGC
  # AGTAAATGCA GATGGTAATA TGCCTTATGA TATTTGTGAT GATGAACAAA CCCTTGACCT TATTGAATCT GAAATGGCTG
  # CTAGAGGAAT TACACAAGAA ATGATTGATG AAAGAAGACA ACAACCAGAA AGGGAAATGT TAAATGATAT GAAAATTTTA
  # CATCAAAGAG GATTACCTTT AGATCAAAGA AATTCTGTTG ATAAATCTAC TTTTGTAAGT TTTTCTGGAG AAAGGGAAAT
  # TTATGTAAAG ATTATTATGA AAGGATTATT ACAGTTTTAT TCCTTTTTAG TTACATATAG CAGCAGCTAA TGGTTATTAT
  # GATGTTGCTG CTTTCCTTCT TCGTTGTAAT GTTTCTCCAG CATTGAGAGA TATAGATTTG TGGCAACCAA TTCATGCAGC
  # TGCTTCTTGG AATCAACCAG ACTTAATCGA GCTTTTATGC GAATATGGGG CTGATATAAA TGCAAAAACT GGAGCTGGGG
  # AAAGCCCTTT AGGTTTATTT TATTGAATCT TATAATTTAT AAATATTTGC TATTAAGTAT GAGGGGAGAG GAACTAACAA
  # TAAGGAATTA AATTTCTCAA TATCAGGATT TTTCGGTTCA CACCCATTTT CTTAAGACCT TTAATTTTTC TCAAAATATG
  # TATGTGACCA CGTCGGGAGG CTTTTTTATT TTTACATGGC TATTTTAAGA AAGGCTAGAA TTTTGACATA CTTTTAACTT
  # ATCGCCTTCC TAACTATTTT CTGTCTATAT ATTTTTTTAA ATTAAGAATT AACTGAAGAT GAACCAACCC AACAAGTAAT
  # TAGAACAATC GCTCAGACAG AAGCAAGGAG ACGGCGTGGT CCAGGTGGTG GTTACTTTGG TGTTCGTGAT TCTCGACGAC
  # AAAGCCGAAA GTAATTTTAA ATTTATATTT TCTTTTCATC TTTTTATCTA GAAGAAAAAA GTTTGAATCT CCTCAACAAC
  # CACCTTCAAC ATTAGAAAAT CCTTTCTCAG CTAGAGGTGC AATTAGACGA CAATCATTGC GAGATCGTAG TGGAATGTCA
  # TTAGCTCGTT TGGAAGCACA AAGAGAGGGT TCTGACCTTA TTAGAAGTTA TAATAGTAAA GAAGACCTTT CTTCTAATAC
  # AGCGGTTTGT TTTTTAAAAT TGTAATTTTT TCTTAATTTT TAGGATGATT CTTTAAATGT TGGAAGTTCT TCATATCTCA
  # ACAATCCAAC AGCCTCGGCT AGTGCTTCCT CTTCAGCATT ACACGGAACT CCACATCAAC AACAACGTCG TGAATCTCCA
  # CCTAAACGTG CATTAATGGC TAGAAGTGCT TCTCATCAAA AACAAAAACA ACAAATGTCT CCAGATGAAT GGCTGAAAAA
  # ATTAGAAGCA GATTCTGCAG GTTTTCGAGA TAATGATGGA GAAGATGGTG AATTACAATC TGAACTTAAA GGAGGACAAA
  # GAATGAAGAG TGGTGGTGGT GGAGGAGCGA GAGGTCAGCA AGGTGAATTA AAATATTTTT TTTGAATTTT ATATTTATTT
  # TTCGTTTAAT AGAAATGAAT GGTGGTCCAA CAGCAACATT TGGTGGAGCT TCAAAACAAC AATTAGCAAT GGGCTCTGGA
  # CCCAATAGAC GGCGCAAACA AGGATGTTGC TCTGTTTTGT GA
  it "should assemble a reverse CDS in MhA1_Contig1133.frz3.gene11" do
    recs = @cdslist['cds:MhA1_Contig1133.frz3.gene11']
    recs.size.should == 8
    component = @componentlist['cds:MhA1_Contig1133.frz3.gene11']
    # 193 bp from MhA1_Contig1133:27,981..28,173
    # >MhA1_Contig1133:27981..28173
    # cgctgtattagaagaaaggtcttctttactattataacttctaataaggtcagaaccctc
    # tctttgtgcttccaaacgagctaatgacattccactacgatctcgcaatgattgtcgtct
    # aattgcacctctagctgagaaaggattttctaatgttgaaggtggttgttgaggagattc
    # aaacttttttctt
    cds1 = recs[5]
    cds1.start.should == 27981
    cds1.frame.should == 1
    cds1.strand.should == '-'
    seq = @gff.assemble(@contigsequence,component.start,[cds1], :codonize=>false)
    seq.should == "CGCTGTATTAGAAGAAAGGTCTTCTTTACTATTATAACTTCTAATAAGGTCAGAACCCTCTCTTTGTGCTTCCAAACGAGCTAATGACATTCCACTACGATCTCGCAATGATTGTCGTCTAATTGCACCTCTAGCTGAGAAAGGATTTTCTAATGTTGAAGGTGGTTGTTGAGGAGATTCAAACTTTTTTCTT"
    seq = @gff.assemble(@contigsequence,component.start,[cds1], :codonize=>true)
    seq.should == "CGCTGTATTAGAAGAAAGGTCTTCTTTACTATTATAACTTCTAATAAGGTCAGAACCCTCTCTTTGTGCTTCCAAACGAGCTAATGACATTCCACTACGATCTCGCAATGATTGTCGTCTAATTGCACCTCTAGCTGAGAAAGGATTTTCTAATGTTGAAGGTGGTTGTTGAGGAGATTCAAACTTTTTTCTT"
    # cds1.frame = 1
    aaseq = @gff.assembleAA(@contigsequence,component.start,[cds1])
    # note it should handle the frame shift and direction!
    # >EMBOSS_001_4
    # RKKFESPQQPPSTLENPFSARGAIRRQSLRDRSGMSLARLEAQREGSDLIRSYNSKEDLSSNTA
    aaseq.should == "RKKFESPQQPPSTLENPFSARGAIRRQSLRDRSGMSLARLEAQREGSDLIRSYNSKEDLSSNTA"
  end
  it "should take the 6th CDS in MhA1_Contig1133.frz3.gene11 (which is 3rd on DNA)" do
    # >MhA1_Contig1133:27981..28173
    # cgctgtattagaagaaaggtcttctttactattataacttctaataaggtcagaaccctc
    # tctttgtgcttccaaacgagctaatgacattccactacgatctcgcaatgattgtcgtct
    # aattgcacctctagctgagaaaggattttctaatgttgaaggtggttgttgaggagattc
    # aaacttttttctt
    recs = @cdslist['cds:MhA1_Contig1133.frz3.gene11']
    component = @componentlist['cds:MhA1_Contig1133.frz3.gene11']
    cds2 = recs[5].clone
    p cds2
    cds2.start.should == 27981
    cds2.frame.should == 1
    cds2.strand.should == '-'
    seq = @gff.assemble(@contigsequence,component.start,[cds2], :codonize=>false)
    seq.should == "CGCTGTATTAGAAGAAAGGTCTTCTTTACTATTATAACTTCTAATAAGGTCAGAACCCTCTCTTTGTGCTTCCAAACGAGCTAATGACATTCCACTACGATCTCGCAATGATTGTCGTCTAATTGCACCTCTAGCTGAGAAAGGATTTTCTAATGTTGAAGGTGGTTGTTGAGGAGATTCAAACTTTTTTCTT"
    aaseq = @gff.assembleAA(@contigsequence,component.start,[cds2])
    # note it should handle the frame shift and direction!
    # >27981..28173_4 RKKFESPQQPPSTLENPFSARGAIRRQSLRDRSGMSLARLEAQREGSDLIRSYNSKEDLSSNTA
    aaseq.should == "RKKFESPQQPPSTLENPFSARGAIRRQSLRDRSGMSLARLEAQREGSDLIRSYNSKEDLSSNTA"
  end
  it "should assemble the 2nd reverse CDS in MhA1_Contig1133.frz3.gene11" do
    recs = @cdslist['cds:MhA1_Contig1133.frz3.gene11']
    component = @componentlist['cds:MhA1_Contig1133.frz3.gene11']
    cds2 = recs[2]
    p cds2
    cds2.frame.should == 1
    cds2.strand.should == '-'
    seq = @gff.assemble(@contigsequence,component.start,[cds2], :codonize=>false)
    seq.should == "CGTTGCGTCTAAATGATCGTCATTTACGTCTACCATTATACGGAATACTATAAACACTACTACTTGTTTGGGAACTGGAATAACTTAGACTTTACCGACGATCTCCTTAATGTGTTCTTTACTAACTACTTTCTTCTGTTGTTGGTCTTTCCCTTTACAATTTACTATACTTTTAAAATGTAGTTTCTCCTAATGGAAATCTAGTTTCTTTAAGACAACTATTTAGATGAAAACATTCAAAAAGACCTCTTTCCCTTTAAATA"
    seq = @gff.assemble(@contigsequence,component.start,[cds2], :codonize=>true)
    seq.should == "GTTGCGTCTAAATGATCGTCATTTACGTCTACCATTATACGGAATACTATAAACACTACTACTTGTTTGGGAACTGGAATAACTTAGACTTTACCGACGATCTCCTTAATGTGTTCTTTACTAACTACTTTCTTCTGTTGTTGGTCTTTCCCTTTACAATTTACTATACTTTTAAAATGTAGTTTCTCCTAATGGAAATCTAGTTTCTTTAAGACAACTATTTAGATGAAAACATTCAAAAAGACCTCTTTCCCTTTAAAT"
    # cds1.frame = 1
    aaseq = @gff.assembleAA(@contigsequence,component.start,[cds2])
    # note it should handle the frame shift and direction!
    aaseq.should == "HGNCYYFLNKQKTNPCLEHNYGFCDLKPVVTCVDVQHEYNYNKLTNNV"
  end
  it "should assemble CDSs, correcting for CODON size"
end

