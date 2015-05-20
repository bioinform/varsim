package com.bina.varsim.fastqLiftover.readers;

import java.util.ArrayList;
import java.io.IOException;
import java.io.BufferedReader;

import com.bina.varsim.fastqLiftover.types.GenomeLocation;
import com.bina.varsim.fastqLiftover.types.MafRecord;
import com.bina.varsim.fastqLiftover.types.SimulatedRead;
import com.bina.varsim.fastqLiftover.types.SimulatedReadPair;
import htsjdk.samtools.fastq.FastqReader;
import htsjdk.samtools.fastq.FastqRecord;

public class PBSIMFastqReader implements PairedFastqReader {

    private FastqReader fastq_;
    private MafReader maf_;
    private Idx2Chr idx2chr_; // base-1 indices to reference name, according to PBSim's enumeration

    public PBSIMFastqReader(final BufferedReader brRef, final BufferedReader brMAF, final BufferedReader brFastq, final boolean forceFiveBaseEncoding) throws IOException {
        maf_ = new MafReader(brMAF);
        fastq_ = new FastqReader(brFastq);
        idx2chr_ = new Idx2Chr(brRef);
    }

    /**
     * get the next entry of PBSim output
     *
     * @return null if no more entries, otherwise a new instance of SimulatedReadPair (read1 is valid and read2 is null)
     */
    @Override
    public SimulatedReadPair getNextReadPair() throws IOException {
        if (fastq_.hasNext()) {
            if ( !maf_.hasNext() ) throw new RuntimeException(); //should formally use zipped iterator

            FastqRecord fastq_entry = fastq_.next();
            MafRecord maf_entry = maf_.next();

            SimulatedRead read = new SimulatedRead();
            read.fragment = 1; //always read-1 since there is no pair-ended-ness
            read.setReadId(fastq_entry.getReadHeader());

            read.sequence = fastq_entry.getReadString();
            read.quality = fastq_entry.getBaseQualityString();

            if (maf_entry.size() != 2) throw new RuntimeException("unexpected MAF data");
            if ( !maf_entry.get(0).src.equals("ref")) throw new RuntimeException("unexpected MAF data");
            if( !maf_entry.get(1).src.equals(read.getReadId())) throw new RuntimeException("unmatched read names");

            //read name is S%d_%d, where the first integer is CHR index and second integer is read number
            final String[] tags = read.getReadId().substring(1).split("_");
            if (tags.length != 2) throw new RuntimeException("unexpected MAF data");

            final GenomeLocation loc = new GenomeLocation(
                    idx2chr_.get(Integer.parseInt(tags[0])),
                    maf_entry.get(0).start0+1, // 0-base to 1-base conversion
                    maf_entry.get(0).strand?0:1); // MAF's "+" maps to 0, "-" to 1

            read.locs1.add(loc);
            read.origLocs1.add(loc);
            read.alignedBases1 = maf_entry.get(1).size;

            return new SimulatedReadPair(read);
        } else {
            return null;
        }
    }

    static private class Idx2Chr {
        private final ArrayList<String> data_;
        public Idx2Chr(final BufferedReader br) throws IOException{
            ArrayList<String> tmp = new ArrayList<String>(5);
            tmp.add(null); // index-0 should not be used, make it break
            for(String buffer = br.readLine() ; buffer != null ; buffer = br.readLine()) {
                tmp.add(buffer.substring(1));
            }
            data_ = tmp;
        }
        String get(int index) {
            return data_.get(index);
        }
    }
}