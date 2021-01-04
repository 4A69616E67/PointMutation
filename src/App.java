import com.github.SnowFlakes.File.FastQFile.*;
import com.github.SnowFlakes.tool.Tools;
import com.github.SnowFlakes.unit.Opts;
import com.github.SnowFlakes.unit.Parameter;
import org.apache.commons.cli.*;
import org.biojava.nbio.alignment.Alignments;
import org.biojava.nbio.alignment.SimpleGapPenalty;
import org.biojava.nbio.core.alignment.matrices.SubstitutionMatrixHelper;
import org.biojava.nbio.core.alignment.template.SequencePair;
import org.biojava.nbio.core.alignment.template.SubstitutionMatrix;
import org.biojava.nbio.core.exceptions.CompoundNotFoundException;
import org.biojava.nbio.core.sequence.DNASequence;
import org.biojava.nbio.core.sequence.compound.AmbiguityDNACompoundSet;
import org.biojava.nbio.core.sequence.compound.NucleotideCompound;


import java.io.*;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;

public class App {
    private final FastqFile Fastq_R1;
    private final FastqFile Fastq_R2;
    private DNASequence Ref_Seq;
    private static final SubstitutionMatrix<NucleotideCompound> Matrix = SubstitutionMatrixHelper.getNuc4_4();
    public int threads = 1;
    public float threshold = 0.99f;
    public String Prefix = "point_mutation";
    public File OutDir = new File("./");
    private final HashMap<String, String[]> BaseMerge = new HashMap<>();

    public App(FastqFile r1, FastqFile r2, String s) {
        Fastq_R1 = r1;
        Fastq_R2 = r2;
        try {
            Ref_Seq = new DNASequence(s, AmbiguityDNACompoundSet.getDNACompoundSet());
        } catch (CompoundNotFoundException e) {
            e.printStackTrace();
        }
        BaseMerge.put("A", new String[]{"A"});
        BaseMerge.put("T", new String[]{"T"});
        BaseMerge.put("C", new String[]{"C"});
        BaseMerge.put("G", new String[]{"G"});
        BaseMerge.put("W", new String[]{"A", "T"});
        BaseMerge.put("M", new String[]{"A", "C"});
        BaseMerge.put("R", new String[]{"A", "G"});
        BaseMerge.put("Y", new String[]{"T", "C"});
        BaseMerge.put("K", new String[]{"T", "G"});
        BaseMerge.put("S", new String[]{"C", "G"});
        BaseMerge.put("N", new String[]{"A", "T", "C", "G"});
    }

    public boolean FileCheck() {
        if ((Fastq_R1 == null || !Fastq_R1.isFile()) && (Fastq_R2 == null || !Fastq_R2.isFile())) {
            System.err.println("error! Incorrect input file:" + Fastq_R1 + "\t" + Fastq_R2);
            return false;
        }
        if (threads <= 0) {
            System.err.println("Error! threads must be Positive");
            return false;
        }
        return true;
    }


    public static void main(String[] args) throws Exception {
        Options Argument = new Options();
        Argument.addOption(Option.builder("1").hasArg().argName("file").desc("fastq R1 file").build());
        Argument.addOption(Option.builder("2").hasArg().argName("file").desc("fastq R2 file").build());
        Argument.addOption(Option.builder("r").hasArg().argName("string").desc("reference seq").required().build());
        Argument.addOption(Option.builder("t").hasArg().argName("int").desc("threads").build());
        Argument.addOption(Option.builder("p").hasArg().argName("string").desc("out prefix").build());
        Argument.addOption(Option.builder("o").hasArg().argName("dir").desc("out dir").build());
        Argument.addOption(Option.builder("f").hasArg().argName("float").desc("threshold, represent how many data will be remained (0< and <1, default 0.99)").build());
        if (args.length <= 0) {
            new HelpFormatter().printHelp("java -jar " + Opts.JarFile.getName(), Argument, true);
            System.exit(1);
        }
        CommandLine comline = new DefaultParser().parse(Argument, args);
        FastqFile f1 = null, f2 = null;
        if (comline.hasOption("1")) {
            f1 = new FastqFile(Parameter.GetStringOpt(comline, "1", ""));
        }
        if (comline.hasOption("2")) {
            f2 = new FastqFile(Parameter.GetStringOpt(comline, "2", ""));
        }
        App app = new App(f1, f2, Parameter.GetStringOpt(comline, "r", null));
        app.threads = Parameter.GetIntOpt(comline, "t", 1);
        app.threshold = Parameter.GetFloatOpt(comline, "f", 0.99f);
        app.OutDir = Parameter.GetFileOpt(comline, "o", new File("./"));
        app.Prefix = Parameter.GetStringOpt(comline, "p", "point_mutation");
        app.run();
    }

    public void run() throws IOException {
        if (!FileCheck()) {
            System.exit(1);
        }
        HashMap<String, int[]> CountHash;
        if (Fastq_R1.isFile() && Fastq_R2.isFile()) {
            CountHash = MultiAlign(Fastq_R1, Fastq_R2, threads);
        } else if (Fastq_R1.isFile()) {
            CountHash = SingleAlign(Fastq_R1, threads);
        } else {
            CountHash = SingleAlign(Fastq_R2, threads);
        }
        //==============================================================================================================
        ArrayList<SeqCount> CountList = Filter(CountHash, threshold);
        File OutFile = new File(OutDir + "/" + Prefix + ".SeqCount.txt");
        BufferedWriter writer = new BufferedWriter(new FileWriter(OutFile));
        writer.write(">target\t" + new DecimalFormat(",###").format(CountList.get(0).TotalCount) + "/100%\n");
        writer.write(CountList.get(0).Seq + "\n");
        for (int i = 1; i < CountList.size(); i++) {
            writer.write(">Seq" + i + "\t" + new DecimalFormat(",###").format(CountList.get(i).TotalCount) + "/" + new DecimalFormat("#.######%").format((float) CountList.get(i).TotalCount / CountList.get(0).TotalCount) + "\n");
            writer.write(CountList.get(i).Seq + "\n");
        }
        writer.close();
        SeqCount ref = CountList.remove(0);
        //--------------------------------------------------------------------------------------------------------------
        String[] BaseList = new String[]{"A", "T", "C", "G"};
        ArrayList<SeqCount> CountStat = Stat(CountList, ref.Seq, BaseList);
        OutFile = new File(OutDir + "/" + Prefix + ".MutationStat.txt");
        writer = new BufferedWriter(new FileWriter(OutFile));
        writer.write(Prefix);
        for (SeqCount seqCount : CountStat) {
            writer.write("\t" + seqCount.Seq);
        }
        writer.write("\n");
        for (int i = 0; i < BaseList.length; i++) {
            writer.write(BaseList[i]);
            for (SeqCount seqCount : CountStat) {
                writer.write("\t" + new DecimalFormat("00.00").format((float) seqCount.Counts[i] / ref.TotalCount * 100));
            }
            writer.write("\n");
        }
        writer.write("Mutation");
        for (SeqCount seqCount : CountStat) {
            writer.write("\t" + new DecimalFormat("00.00").format((float) seqCount.TotalCount / ref.TotalCount * 100));
        }
        writer.write("\n");
        writer.close();
    }

    private HashMap<String, int[]> SingleAlign(FastqFile fastqFile, int threads) throws IOException {
        HashMap<String, int[]> SeqCount = new HashMap<>();
        FastqFile inFile = new FastqFile(fastqFile.getPath());
        inFile.ReadOpen();
        Thread[] t = new Thread[threads];
        for (int i = 0; i < t.length; i++) {
            t[i] = new Thread(() -> {
                try {
                    FastqItem item;
                    DNASequence target;
                    while ((item = fastqFile.ReadItem()) != null) {
                        target = new DNASequence(item.Sequence, AmbiguityDNACompoundSet.getDNACompoundSet());
                        SequencePair<DNASequence, NucleotideCompound> psa1 = Alignments.getPairwiseAlignment(Ref_Seq, target, Alignments.PairwiseSequenceAlignerType.LOCAL, new SimpleGapPenalty(10, 2), Matrix);
                        SequencePair<DNASequence, NucleotideCompound> psa2 = Alignments.getPairwiseAlignment(Ref_Seq, new DNASequence(target.getReverseComplement().getSequenceAsString(), AmbiguityDNACompoundSet.getDNACompoundSet()), Alignments.PairwiseSequenceAlignerType.LOCAL, new SimpleGapPenalty(10, 2), Matrix);
                        SequencePair<DNASequence, NucleotideCompound> psa = psa1.getNumSimilars() >= psa2.getNumSimilars() ? psa1 : psa2;
                        String s = psa.getAlignedSequence(2).toString().replaceAll("-", "");
                        synchronized ("lock1") {
                            if (!SeqCount.containsKey(s)) {
                                SeqCount.put(s, new int[]{0});
                            }
                            SeqCount.get(s)[0]++;
                        }
                    }
                } catch (Exception e) {
                    e.printStackTrace();
                }
            });
            t[i].start();
        }
        Tools.ThreadsWait(t);
        inFile.ReadClose();
        return SeqCount;
    }

    private HashMap<String, int[]> MultiAlign(FastqFile fastq_R1, FastqFile fastq_R2, int threads) throws IOException {
        HashMap<String, int[]> SeqCount = new HashMap<>();
        FastqFile r1 = new FastqFile(fastq_R1.getPath());
        FastqFile r2 = new FastqFile(fastq_R2.getPath());
        r1.ReadOpen();
        r2.ReadOpen();
        Thread[] t = new Thread[threads];
        for (int i = 0; i < t.length; i++) {
            t[i] = new Thread(() -> {
                try {
                    FastqItem item1;
                    FastqItem item2;
                    synchronized ("lock2") {
                        item1 = r1.ReadItem();
                        item2 = r2.ReadItem();
                    }
                    DNASequence target = null;
                    while (item1 != null && item2 != null) {
                        if (item1.Title.split("\\s+")[0].equals(item2.Title.split("\\s+")[0])) {
                            DNASequence DS1 = new DNASequence(item1.Sequence, AmbiguityDNACompoundSet.getDNACompoundSet());
                            DNASequence DS2 = new DNASequence(new DNASequence(item2.Sequence).getReverseComplement().getSequenceAsString(), AmbiguityDNACompoundSet.getDNACompoundSet());
                            SequencePair<DNASequence, NucleotideCompound> psa = Alignments.getPairwiseAlignment(DS1, DS2, Alignments.PairwiseSequenceAlignerType.GLOBAL, new SimpleGapPenalty(10, 0), Matrix);
                            target = Merge(psa);
                        } else {
                            System.err.println("Different title in R1 and R2 file:\tR1:" + item1.Title + "\tR2:" + item2.Title);
                            System.exit(1);
                        }
                        if (target != null) {
                            SequencePair<DNASequence, NucleotideCompound> psa1 = Alignments.getPairwiseAlignment(Ref_Seq, target, Alignments.PairwiseSequenceAlignerType.LOCAL, new SimpleGapPenalty(10, 2), Matrix);
                            SequencePair<DNASequence, NucleotideCompound> psa2 = Alignments.getPairwiseAlignment(Ref_Seq, new DNASequence(target.getReverseComplement().getSequenceAsString(), AmbiguityDNACompoundSet.getDNACompoundSet()), Alignments.PairwiseSequenceAlignerType.LOCAL, new SimpleGapPenalty(10, 2), Matrix);
                            SequencePair<DNASequence, NucleotideCompound> psa = psa1.getNumSimilars() >= psa2.getNumSimilars() ? psa1 : psa2;
                            String s = psa.getAlignedSequence(2).toString().replaceAll("-", "");
                            synchronized ("lock1") {
                                if (!SeqCount.containsKey(s)) {
                                    SeqCount.put(s, new int[]{0});
                                }
                                SeqCount.get(s)[0]++;
                            }
                        }
                        synchronized ("lock2") {
                            item1 = r1.ReadItem();
                            item2 = r2.ReadItem();
                        }
                    }
                } catch (Exception e) {
                    e.printStackTrace();
                }
            });
            t[i].start();
        }
        Tools.ThreadsWait(t);
        r1.ReadClose();
        r2.ReadClose();
        return SeqCount;
    }

    private ArrayList<SeqCount> Filter(HashMap<String, int[]> map, float threshold) {
        ArrayList<SeqCount> CountList = new ArrayList<>();
        int TotalNum = 0;
        for (String key : map.keySet()) {
            CountList.add(new SeqCount(key, map.get(key)[0]));
            TotalNum += map.get(key)[0];
        }
        float ThresholdVal = TotalNum * threshold;
        CountList.sort(new SeqCount.CountComparator());
        Collections.reverse(CountList);
        int i;
        TotalNum = 0;
        for (i = 0; i < CountList.size(); i++) {
            TotalNum += CountList.get(i).TotalCount;
            if (TotalNum >= ThresholdVal) {
                break;
            }
        }
        if (CountList.size() > i + 1) {
            CountList.subList(i + 1, CountList.size()).clear();
        }
        CountList.add(0, new SeqCount(this.Ref_Seq.getSequenceAsString(), TotalNum));
        return CountList;
    }

    private ArrayList<SeqCount> Stat(ArrayList<SeqCount> arrayList, String ref, String[] bases) {
        ArrayList<SeqCount> CountStat = new ArrayList<>();
        for (String c : ref.split("")) {
            SeqCount item = new SeqCount(c, 0);
            item.Counts = new int[bases.length];
            CountStat.add(item);
        }
        CountStat.add(new SeqCount("Insert", new int[bases.length], 0));
        CountStat.add(new SeqCount("Delete", new int[bases.length], 0));
        for (SeqCount item : arrayList) {
            if (item.Seq.length() > ref.length()) {
                CountStat.get(ref.length()).TotalCount += item.TotalCount;
            } else if (item.Seq.length() < ref.length()) {
                CountStat.get(ref.length() + 1).TotalCount += item.TotalCount;
            } else {
                String[] item_seq = item.Seq.split("");
                for (int j = 0; j < item.Seq.length(); j++) {
                    for (int k = 0; k < bases.length; k++) {
                        if (item_seq[j].compareToIgnoreCase(bases[k]) == 0) {
                            CountStat.get(j).Counts[k] += item.TotalCount;
                            break;
                        }
                    }
                }
            }
        }
        for (int i = 0; i < CountStat.size() - 2; i++) {
            String[] ref_bases = this.BaseMerge.get(CountStat.get(i).Seq.toUpperCase());
            for (int j = 0; j < bases.length; j++) {
                boolean flag = true;
                for (String ref_base : ref_bases) {
                    if (bases[j].compareToIgnoreCase(ref_base) == 0) {
                        flag = false;
                        break;
                    }
                }
                if (flag) {
                    CountStat.get(i).TotalCount += CountStat.get(i).Counts[j];
                }
            }

        }
        return CountStat;
    }


    /**
     * @param psa align result with R1 seq and R2 seq
     * @return merged string or null if can't merge
     */
    private DNASequence Merge(SequencePair<DNASequence, NucleotideCompound> psa) throws CompoundNotFoundException {
        char[] res1 = psa.getAlignedSequence(1).toString().toCharArray();
        char[] res2 = psa.getAlignedSequence(2).toString().toCharArray();
        char[] res = new char[res1.length];
        for (int i = 0; i < res1.length; i++) {
            if (res1[i] == '-' || (res1[i] == 'N' && res2[i] != '-') || res1[i] == res2[i]) {
                res[i] = res2[i];
            } else if (res2[i] == '-' || res2[i] == 'N') {
                res[i] = res1[i];
            } else {
                return null;
            }
        }
        return new DNASequence(new String(res), AmbiguityDNACompoundSet.getDNACompoundSet());
    }

//
}

class SeqCount {
    public String Seq;
    public int[] Counts;
    public int TotalCount;

    public SeqCount(String s, int c) {
        Seq = s;
        TotalCount = c;
    }

    public SeqCount(String s, int[] cs, int c) {
        this(s, c);
        Counts = cs;
    }

    static class CountComparator implements Comparator<SeqCount> {

        @Override
        public int compare(SeqCount o1, SeqCount o2) {
            return o1.TotalCount - o2.TotalCount;
        }
    }
}


