package com.bina.varsim.fastqLiftover.readers;

import com.bina.varsim.fastqLiftover.types.GenomeLocation;
import com.bina.varsim.fastqLiftover.types.MafRecord;
import com.bina.varsim.fastqLiftover.types.SimulatedRead;
import com.bina.varsim.fastqLiftover.types.SimulatedReadPair;
import com.bina.varsim.types.ChrString;
import htsjdk.samtools.fastq.FastqReader;
import htsjdk.samtools.fastq.FastqRecord;

import java.io.BufferedReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

public class PBSIMFastqReader implements PairedFastqReader {

    private FastqReader fastq;
    private MafReader maf;
    private Idx2Chr idx2Chr; // base-1 indices to reference name, according to PBSim's enumeration

    public PBSIMFastqReader(final BufferedReader brRef, final BufferedReader brMAF, final BufferedReader brFastq, final boolean forceFiveBaseEncoding) throws IOException {
        maf = new MafReader(brMAF);
        fastq = new FastqReader(brFastq);
        idx2Chr = new Idx2Chr(brRef);
    }

    /**
     * get the next entry of PBSim output
     *
     * @return null if no more entries, otherwise a new instance of SimulatedReadPair (read1 is valid and read2 is null)
     */
    @Override
    public SimulatedReadPair getNextReadPair() throws IOException {
        if (fastq.hasNext()) {
            if (!maf.hasNext()) throw new RuntimeException(); //should formally use zipped iterator

            FastqRecord fastqEntry = fastq.next();
            MafRecord mafEntry = maf.next();

            SimulatedRead read = new SimulatedRead();
            read.fragment = 1; //always read-1 since there is no pair-ended-ness
            read.setReadId(fastqEntry.getReadHeader());

            read.sequence = fastqEntry.getReadString();
            read.quality = fastqEntry.getBaseQualityString();

            if (mafEntry.size() != 2) throw new RuntimeException("unexpected MAF data");
            if (!mafEntry.get(0).src.equals("ref")) throw new RuntimeException("unexpected MAF data");
            if (!mafEntry.get(1).src.equals(read.getReadId())) throw new RuntimeException("unmatched read names");

            //read name is S%d_%d, where the first integer is CHR index and second integer is read number
            final String[] tags = read.getReadId().substring(1).split("_");
            if (tags.length != 2) throw new RuntimeException("unexpected MAF data");

            final GenomeLocation loc = new GenomeLocation(
                    idx2Chr.get(Integer.parseInt(tags[0])),
                    mafEntry.get(0).start0 + 1, // 0-base to 1-base conversion
                    mafEntry.get(0).strand ? 0 : 1); // MAF's "+" maps to 0, "-" to 1

            read.locs1.add(loc);
            read.origLocs1.add(loc);
            read.alignedBases1 = mafEntry.get(1).size;

            return new SimulatedReadPair(read);
        } else {
            return null;
        }
    }

    static private class Idx2Chr {
        private final List<ChrString> data;

        public Idx2Chr(final BufferedReader br) throws IOException {
            List<ChrString> tmp = new ArrayList<>(5);
            tmp.add(null); // index-0 should not be used, make it break
            for (String buffer = br.readLine(); buffer != null; buffer = br.readLine()) {
                tmp.add(new ChrString(buffer.substring(1)));
            }
            data = tmp;
        }

        ChrString get(int index) {
            return data.get(index);
        }
    }
}