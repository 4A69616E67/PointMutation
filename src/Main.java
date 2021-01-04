import com.github.SnowFlakes.File.BedFile.BedFile;
import com.github.SnowFlakes.File.BedFile.BedItem;
import com.github.SnowFlakes.File.FastQFile.FastqFile;
import com.github.SnowFlakes.Software.Bwa;
import com.github.SnowFlakes.unit.Opts;
import com.github.SnowFlakes.unit.Parameter;
import org.apache.commons.cli.*;

import java.io.*;
import java.util.ArrayList;

public class Main {
    /**
     * Format
     * chr  start   end 0   +   name    seq
     */
    private BedFile LocationFile;
    private FastqFile Fastq1, Fastq2;
    private File IndexGenome;
    private Bwa bwa = new Bwa("bwa");
    public int Threads = 1;

    public Main(FastqFile fq1, FastqFile fq2, File index, BedFile locationFile) {
        Fastq1 = fq1;
        Fastq2 = fq2;
        IndexGenome = index;
        LocationFile = locationFile;
    }

    public boolean FileCheck() {
        if (!Fastq1.exists() || Fastq1.isDirectory()) {
            System.err.println("Error! Incorrect R1 file:\t" + Fastq1);
            return false;
        }
        if (!Fastq2.exists() || Fastq2.isDirectory()) {
            System.err.println("Waring! Incorrect or no R2 file:\t" + Fastq2);
        }
        bwa.IndexPrefix = IndexGenome;
        if (!bwa.IndexCheck()) {
            System.err.println("Error! Incorrect bwa index prefix:\t" + IndexGenome);
            return false;
        }
        if (!LocationFile.exists()) {
            System.err.println("Error! Incorrect location file:\t" + LocationFile);
            return false;
        }
        return true;
    }

    public static void main(String[] args) throws IOException {
        Options Argument = new Options();
        Argument.addOption(Option.builder("1").hasArg().argName("file").desc("fastq R1 file").required().build());
        Argument.addOption(Option.builder("2").hasArg().argName("file").desc("fastq R2 file").build());
        Argument.addOption(Option.builder("l").hasArg().argName("file").desc("location file (bed format 'chr start end 0 + name seq')").required().build());
        Argument.addOption(Option.builder("i").hasArg().argName("file").desc("bwa index prefix").required().build());
        Argument.addOption(Option.builder("t").hasArg().argName("int").desc("threads").build());
        if (args.length <= 0) {
            new HelpFormatter().printHelp("java -jar " + Opts.JarFile.getName(), Argument, true);
            System.exit(1);
        }
        CommandLine comline = null;
        try {
            comline = new DefaultParser().parse(Argument, args);
        } catch (ParseException e) {
            System.err.println(e.getMessage());
            System.exit(1);
        }
        FastqFile fq1 = new FastqFile(Parameter.GetStringOpt(comline, "1", ""));
        FastqFile fq2 = new FastqFile(Parameter.GetStringOpt(comline, "2", ""));
        File indexFile = Parameter.GetFileOpt(comline, "i", new File(""));
        BedFile locationFile = new BedFile(Parameter.GetStringOpt(comline, "l", ""));
        Main main = new Main(fq1, fq2, indexFile, locationFile);
        main.Threads = Parameter.GetIntOpt(comline, "t", 1);
        main.run();
    }

    public void run() throws IOException {
        if (!FileCheck()){
            System.exit(1);
        }
        ArrayList<BedItem> locationList = ReadLocation(LocationFile);

    }


    private ArrayList<BedItem> ReadLocation(File locationFile) throws IOException {
        ArrayList<BedItem> list = new ArrayList<>();
        BedFile inFile = new BedFile(locationFile.getPath());
        inFile.ReadOpen();
        BedItem item;
        while ((item = inFile.ReadItem()) != null) {
            list.add(item);
        }
        return list;
    }

}
