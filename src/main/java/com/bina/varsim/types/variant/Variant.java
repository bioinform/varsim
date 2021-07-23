package com.bina.varsim.types.variant;

//--- Java imports ---

import com.bina.intervalTree.SimpleInterval1D;
import com.bina.varsim.types.*;
import com.bina.varsim.types.variant.alt.Alt;
import com.bina.varsim.util.SimpleReference;
import com.bina.varsim.util.StringUtilities;
import org.apache.log4j.Logger;

import java.io.PrintWriter;
import java.io.UnsupportedEncodingException;
import java.util.*;
import java.util.stream.Collectors;

import static java.lang.Character.getType;
import static java.lang.Math.abs;

public class Variant implements Comparable<Variant>{
    private final static Logger log = Logger.getLogger(Variant.class.getName());
    public int splitVariantIndex = 0; // this is hopefully a unique index, for the split variants
    public int wholeVariantIndex = 0; // this is hopefully a unique index, for the whole variants
    // this is the type before the variant was split into canonical ones
    public VariantOverallType originalType = null;
    // use a seed for reproducibility, should be an option or global
    private Random rand = null;
    private int pos = -1, referenceAlleleLength = -1;
    private byte[] ref;
    private ChrString chr;
    private Alt[] alts;
    private byte maternal = 0, paternal = 0; // -1 for not avaliable
    private boolean isPhased = false; // Phasing
    private String filter;
    private String qual;
    private String info;
    private String varId;
    // this is when the reference base is deleted
    // if it is the same as the first alt base
    //so refDeleted.length() <= 1 is always true?
    private String refDeleted;
    private String extraBase = "";
    private String clippedSequence = "";
    private ChrString[] chr2;
    private int[] pos2;
    private int[] end2;
    private Boolean isinv; //is sequence inverted? useful for interspersed dup, translocation dup
    private String traid; //translocation ID
    private volatile int hashCode; //caching hashcode
    private List<Variant> compositions; //when this variant is a composite variant, this field will store constitutional variants

    //additional INFO fields
    private int threePrimeDistance = -1; //3' end distance with a matching variant
    private int fivePrimeDistance = -1;//5' end distance with a matching variant
    private int lengthDifference = -1; //length difference with a matching variant
    private Boolean isLengthImprecise = false; //true if length is imprecise


    public List<Variant> getCompositions() {
        return compositions;
    }

    public void setCompositions(List<Variant> compositions) {
        this.compositions = compositions;
    }



    public Variant(final Random rand) {
        // TODO define some methods to determine if a Variant is uninitialised
        this.rand = rand;
    }

    /**
     * whether have flag ISINV?
     * @return
     */
    public boolean isInversed() {
        return isinv != null && isinv.booleanValue();
    }

    /*
    use Builder pattern to make code more readable and maintainable.
     */
    public static class Builder {
        private Random rand = null;
        private int pos = -1, referenceAlleleLength = -1;
        private byte[] ref;
        private ChrString chr;
        private Alt[] alts;
        private byte maternal = 0, paternal = 0; // -1 for not avaliable
        private boolean isPhased = false; // Phasing
        private String filter;
        private String qual;
        private String info;
        private String varId;
        // this n the reference base is deleted
        // if ite same as the first alt base
        //so refd.length() <= 1 is always true?
        private String refDeleted;
        private String clippedSequence = "";
        private ChrString[] chr2;
        private int[] pos2;
        private int[] end2;
        private Boolean isinv; //is sequence inverted? useful for interspersed dup, translocation dup
        private Boolean isLengthImprecise = false;
        private String traid; //translocation ID
        private List<Variant> compositions;

        public Builder() {}
        public Builder chr(final ChrString chr) {
            this.chr = chr;
            return this;
        }
        public Builder pos(final int pos) {
            this.pos = pos;
            return this;
        }
        public Builder referenceAlleleLength(final int referenceAlleleLength) {
            this.referenceAlleleLength = referenceAlleleLength;
            return this;
        }
        public Builder ref(final byte[] ref) {
            this.ref = ref.clone();
            return this;
        }
        public Builder alts(final Alt[] alts) {
            this.alts = new Alt[alts.length];
            for (int i = 0; i < alts.length; i++) {
                if (alts[i] != null) {
                    this.alts[i] = alts[i].copy(); //deep copy
                } else {
                    this.alts[i] = null;
                }
            }
            return this;
        }
        public Builder phase(final byte[] phase) {
            this.paternal = phase[0];
            this.maternal = phase[1];
            return this;
        }
        public Builder isPhased(final boolean isPhased) {
            this.isPhased = isPhased;
            return this;
        }
        public Builder varId(final String varId) {
            this.varId = varId;
            return this;
        }
        public Builder filter(final String filter) {
            this.filter = filter;
            return this;
        }
        public Builder qual(final String qual) {
            this.qual = qual;
            return this;
        }
        public Builder info(final String info) {
            this.info = info;
            return this;
        }
        public Builder refDeleted(final String refDeleted) {
            this.refDeleted = refDeleted;
            return this;
        }
        public Builder randomNumberGenerator(final Random rand) {
            this.rand = rand;
            return this;
        }
        public Builder chr2(final ChrString[] chr2) {
            this.chr2 = chr2;
            return this;
        }
        public Builder pos2(final int[] pos2) {
            this.pos2 = pos2;
            return this;
        }
        public Builder end2(final int[] end2) {
            this.end2 = end2;
            return this;
        }
        public Builder isinv(final Boolean b) {
            this.isinv = b;
            return this;
        }
        public Builder isLengthImprecise(final Boolean b) {
            this.isLengthImprecise = b;
            return this;
        }
        public Builder traid(final String id) {
            this.traid = id;
            return this;
        }
        public Builder compositions(final List<Variant> compositions) {
            this.compositions = compositions;
            return this;
        }
        public Builder clippedSequence(final String clippedSequence) {
            this.clippedSequence = clippedSequence;
            return this;
        }
        public Variant build() {
            return new Variant(this);
        }
    }
    private Variant(Builder builder) {
        this.rand = builder.rand;
        this.varId = builder.varId;
        this.filter = builder.filter;
        this.qual = builder.qual;
        this.info = builder.info;
        this.chr = builder.chr;
        this.pos = builder.pos;
        this.referenceAlleleLength = builder.referenceAlleleLength;

        this.ref = builder.ref;
        this.refDeleted = builder.refDeleted;
        this.clippedSequence = builder.clippedSequence;
        this.alts = builder.alts;
        this.chr2 = builder.chr2;
        this.pos2 = builder.pos2;
        this.end2 = builder.end2;
        this.paternal = builder.paternal;
        this.maternal = builder.maternal;
        this.isPhased = builder.isPhased;
        this.traid = builder.traid;
        this.isinv = builder.isinv;
        this.isLengthImprecise = builder.isLengthImprecise;
        this.compositions = builder.compositions;
    }

    /*consider making Alt class immutable, all breakend, symbolic allele interpretation,
    alt allele sequence accessing/storing handled internally and automatically
    without explicitly invoking any public methods.

    also make the Variant constructor with many arguments private, so that clients
    will not use it any more. this will force the client to use Builder and then
    we are free to change the constructor.
    */

    public Variant(final Variant var) {
        filter = var.filter;
        qual = var.qual;
        info = var.info;
        varId = var.varId;
        chr = var.chr;
        pos = var.pos;
        referenceAlleleLength = var.referenceAlleleLength;
        ref = var.ref == null ? null : var.ref.clone();
        refDeleted = var.refDeleted;
        clippedSequence = var.clippedSequence;
        if (var.alts != null) {
            alts = new Alt[var.alts.length];
            for (int i = 0; i < var.alts.length; i++) {
                if (var.alts[i] != null) {
                    alts[i] = var.alts[i].copy();
                } else {
                    alts[i] = null;
                }
            }
        }

        this.chr2 = var.chr2;
        this.pos2 = var.pos2;
        this.end2 = var.end2;
        paternal = var.paternal;
        maternal = var.maternal;
        isPhased = var.isPhased;
        rand = var.rand;
        isinv = var.isinv;
        isLengthImprecise = var.isLengthImprecise;
        traid = var.traid;
        compositions = var.getCompositions();
    }

    /**
     * @return Chromosome variant is on
     */
    public ChrString getChr() {
        return chr;
    }

    /**
     * @return Start position of variant
     */
    public int getPos() {
        return pos;
    }

    // return false if fails
    // if return false, nothing is changed

    /**
     * Tries to move the variant to a new novel position. It checks if the new position is valid
     *
     * @param pos new novel position
     * @param ref reference sequence
     * @return return true if it was changed
     */
    public boolean setNovelPosition(final int pos, final SimpleReference ref) {

        // replace ref
        int len = this.ref.length;

        if (len > 0) {
            byte[] temp_ref = ref.byteRange(chr, pos, pos + len);

            for (byte b : temp_ref) {
                if (b == 'N') {
                    // don't allow N's
                    log.warn("N found at " + pos + " to " + (pos + len));
                    return false;
                }
            }

            for (Alt alt : alts) {
                if (alt.getSequence() != null) {
                    // make sure there is no prefix the same
                    for (int i = 0; i < temp_ref.length; i++) {
                        if (i < alt.getSequence().length) {
                            if (temp_ref[i] == alt.getSequence()[i]) {
                                log.warn("Same ref at alt at " + pos + " to " + (pos + len));
                                return false;
                            } else {
                                break;
                            }
                        }
                    }

                }
            }

            this.ref = temp_ref;
        }

        // replace ref_deleted
        len = refDeleted.length();
        if (len > 0) {
            try {
                byte[] deleted_temp = ref.byteRange(chr, pos - len, pos);
                if (deleted_temp == null) {
                    throw new IllegalArgumentException("Error record: " + this);
                }
                refDeleted = new String(deleted_temp, "US-ASCII");
            } catch (UnsupportedEncodingException e) {
                e.printStackTrace();
            }
        }

        // do the replacing
        this.pos = pos;

        return true;
    }

    /**
     * @param id variant id. usually the dbSNP id
     */
    public void setVarID(final String id) {
        varId = id;
    }

    /**
     * return the length of the getReferenceAlleleLength.
     * This is currently simply the length of the reference sequence that will
     * be replaced
     *
     * @return the length of the getReferenceAlleleLength
     */
    public int getReferenceAlleleLength() {
        return referenceAlleleLength;
    }

    /**
     * @return maternal allele index, 0 if reference
     */
    public int maternal() {
        return maternal;
    }

    /**
     * @return paternal allele index, 0 if reference
     */
    public int paternal() {
        return paternal;
    }

    /**
     * @param ind index of allele
     * @return the insertion sequence as a string
     */
    public byte[] insertion(final int ind) {
        return (ind <= 0 || ind > alts.length) ? null : alts[ind - 1].getSequence();
    }

    public ChrString getChr2(final int ind) {
        return (ind <= 0 || ind > alts.length) ? new ChrString("") : this.chr2[ind - 1];
    }

    public int getPos2(final int ind) {
        return (ind <= 0 || ind > alts.length) ? -1 : this.pos2[ind - 1];
    }

    public int getEnd2(final int ind) {
        return (ind <= 0 || ind > alts.length) ? -1 : this.end2[ind - 1];
    }

    public ChrString[] getAllChr2() {
        return this.chr2 == null ? new ChrString[0] : this.chr2;
    }

    public int[] getAllPos2() {
        return this.pos2 == null ? new int[0] : this.pos2;
    }

    public int[] getAllEnd2() {
        return this.end2 == null ? new int[0] : this.end2;
    }

    /**
     * return 1-based end of variant depending on overall type of the variant
     * for insertion, SNP, interspersed duplication, tandem duplication
     * and translocation duplication, return 3' end
     * otherwise, return 3' end plus length
     * @return
     */
    public int getEnd() {
        switch(getType()) {
            case Insertion:
            case SNP:
            case InterDup:
            case TandemDup:
            case TransDup:
                return pos;
            default:
                return pos + maxLen() - 1;
        }
    }

    public String getTraid() {
        return traid;
    }

    /**
     * The length of an alternate allele, this is usually an insertion sequence
     * But in the case of SVs, not necessarily
     *
     * @param ind index of allele
     * @return the length of that allele
     */
    public int getInsertionLength(final int ind) {
        return (ind <= 0 || ind > alts.length) ? 0 : alts[ind - 1].length();
    }

    /** if it is a simple indel, it is just the length
     if it is a complex variant, this is the maximum length of the insertion
     and getReferenceAlleleLength

     @param ind index of allele
     */
    public int maxLen(final int ind) {
        return (ind <= 0 || ind > alts.length) ? 0 : Math.max(referenceAlleleLength, alts[ind - 1].length());
    }

    /**
     * get maximum length of a variant across
     * genotypes and reference and alternative alleles
     * @return
     */
    public int maxLen() {
        if (compositions == null) {
            return Math.max(maxLen(getAllele(0)), maxLen(getAllele(1)));
        } else {
            if (compositions.isEmpty())
                return 0;
            int maxLen = compositions.get(0).maxLen();
            for (Variant c : compositions) {
                maxLen = Math.max(maxLen, c.maxLen());
            }
            return maxLen;
        }
    }

    // this is the minimum length of the variants
    public int minLen() {
        if (compositions == null) {
            return Math.min(maxLen(getAllele(0)), maxLen(getAllele(1)));
        } else {
            if (compositions.isEmpty())
                return 0;
            int minLen = compositions.get(0).minLen();
            for (Variant c : compositions) {
                minLen = Math.min(minLen, c.minLen());
            }
            return minLen;
        }
    }

    /*
    gets the interval enclosing the variant on the reference genome
    */
    public SimpleInterval1D getAlternativeAlleleInterval(final int ind) {
        if (ind == 0 || referenceAlleleLength == 0) {
            return new SimpleInterval1D(pos, pos);
        }

        return new SimpleInterval1D(pos, pos + referenceAlleleLength - 1);
    }

    /*
    gets the interval for the variant, accounting for variant size
    */
    public SimpleInterval1D getVariantInterval(final int ind) {
        try {
            if (ind == 0) {
                return new SimpleInterval1D(pos, pos);
            }

            // TODO hmm unsafe
            if (maxLen(ind) == Integer.MAX_VALUE) {
                return new SimpleInterval1D(pos, pos);
            } else {
                return new SimpleInterval1D(pos, pos + maxLen(ind) - 1);
            }
        } catch (RuntimeException e) {
            log.error("Bad variant interval: " + toString());
            log.error("pos: " + pos);
            log.error("ind: " + ind);
            log.error("maxLen(ind): " + maxLen(ind));
            e.printStackTrace();
            return null;
        }
    }

    /**
     * if the variant is an insertion and ignoreInsertionLength is true
     * then return interval of length 1, otherwise, return normal
     * interval
     *
     * @param ind
     * @return
     */
    public SimpleInterval1D getVariantInterval(final int ind, final boolean ignoreInsertionLength) {
        if (getType(ind) == VariantType.Insertion && ignoreInsertionLength) {
            return new SimpleInterval1D(getPos(), getPos());
        } else {
            return getVariantInterval(ind);
        }
    }

    // union of intervals from the genotypes
    public SimpleInterval1D getGenotypeUnionAlternativeInterval() {
        return getAlternativeAlleleInterval(getGoodPaternal()).union(getAlternativeAlleleInterval(getGoodMaternal()));
    }

    public SimpleInterval1D getGenotypeUnionVariantInterval() {
      if (getGoodPaternal() == -1 && getGoodMaternal() == -1) {
          throw new RuntimeException("Both maternal and paternal genotypes are missing for " + this.toString());
      } else if (getGoodPaternal() == -1) {
          return getVariantInterval(getGoodMaternal());
      } else if (getGoodMaternal() == -1) {
          return getVariantInterval(getGoodPaternal());
      } else {
          return getVariantInterval(getGoodPaternal()).union(getVariantInterval(getGoodMaternal()));
      }
    }

    /**
     * decide if this variant is a small variant or not
     * @param genotype
     * @param cutoff
     * @return
     */
    public boolean isSmallVariant(int genotype, int cutoff, boolean ignoreInsertionLength) {
        VariantType type = this.getType(genotype);
        SimpleInterval1D intervalForCompare = this.getVariantInterval(genotype, ignoreInsertionLength);
        return ((type ==  VariantType.Insertion || type == VariantType.Deletion || type == VariantType.Complex ) &&
                intervalForCompare.right - intervalForCompare.left + 1 < cutoff);
    }

    public Genotypes getGenotypes() {
        return new Genotypes(getGoodPaternal(), getGoodMaternal());
    }

    /*
    * 0 = paternal
    * 1 = maternal
    * otherwise returns -1
     */
    public int getAllele(final int parent) {
        if (parent == 0) {
            return getGoodPaternal();
        } else if (parent == 1) {
            return getGoodMaternal();
        }
        return -1;
    }

    /**
     * find the length common prefix between reference and alternative alleles
     * @return
     */
    public int getMinMatchLength() {
        int[] allele = {getAllele(0), getAllele(1)};
        byte[][] alt = {getAlt(allele[0]).getSequence(), getAlt(allele[1]).getSequence()};
        byte[] ref = getReference();

        int[] matchLength = {0, 0};
        for (int i = 0; i < 2; i++) {
            for (int j = 0; j < Math.min(ref.length, alt[i].length); j++) {
                if (alt[i][j] == ref[j]) {
                    matchLength[i]++;
                } else {
                    break;
                }
            }
        }
        return Math.min(matchLength[0], matchLength[1]);
    }

    public void setAllele(final int parent, final byte allele) {
        if (parent == 0) {
            paternal = allele;
        } else if (parent == 1) {
            maternal = allele;
        }
    }

    // TODO this is wrong, but it only effects the count of variant bases
    public int variantBases() {
        int ret = referenceAlleleLength;
        for (Alt alt : alts) {
            if (referenceAlleleLength != alt.length()) {
                ret += alt.length();
            }
        }
        return ret;
    }

    /**
     * @param ind index of allele (starts at 1, 0 is reference)
     * @return type of allele at index ind
     */
    public VariantType getType(final int ind) {
        if (ind <= 0) {
            return VariantType.Reference;
        }

        if (compositions != null) {
            return VariantType.Composite;
        }

        /*when variant type is explicitly specified in VCF,
        alt allele will be set to the correct type,
        we can be assured that correct variant type is
        returned.
         */
        Alt alt = null;
        try {
            alt = alts[ind - 1];
        } catch (IndexOutOfBoundsException e) {
            throw e;
        }
        if (alt.getSymbolicAllele() != null) {
            Alt.SVType major = alt.getSymbolicAllele().getMajor();
            Alt.SVType.SVSubtype minor = alt.getSymbolicAllele().getMinor();
            if (minor == Alt.SVType.SVSubtype.TRA) {
                return major == Alt.SVType.DUP ? VariantType.Translocation_Duplication : VariantType.Translocation_Deletion;
            } else if (major == Alt.SVType.INS) {
                return VariantType.Insertion;
            } else if (major == Alt.SVType.DEL && (minor == null || minor != Alt.SVType.SVSubtype.TRA)) {
                return VariantType.Deletion;
            } else if (major == Alt.SVType.INV) {
                return VariantType.Inversion;
            } else if (major == Alt.SVType.DUP && (minor == null || minor == Alt.SVType.SVSubtype.TANDEM)) {
                return VariantType.Tandem_Duplication;
            } else if (major == Alt.SVType.DUP && (minor == Alt.SVType.SVSubtype.ISP || (minor == null && pos2 != null))) {
                return VariantType.Interspersed_Duplication;
            } else {
                return VariantType.Complex;
            }
        } else if (alt.getBreakend() != null) {
            return VariantType.Breakend;
        }
        /*
        FlexSeq is the least reliable source of variant type
        because some variants will be transformed into different
        seq types
         */
        FlexSeq.Type type = alts[ind - 1].getSeqType();
        switch (type) {
            case TANDEM_DUP:
                return VariantType.Tandem_Duplication;
            case INS:
                return VariantType.Insertion;
            case INV:
                return VariantType.Inversion;
            case TRA_DEL:
                return VariantType.Translocation_Deletion;
            case TRA_DUP:
                return VariantType.Translocation_Duplication;
            case ISP_DUP:
                return VariantType.Interspersed_Duplication;
            case DEL:
                return VariantType.Deletion;
            default:
                break;
        }

        int insertionLength = getInsertionLength(ind);
        int deletionLength = referenceAlleleLength;
        if (insertionLength == 0 && deletionLength == 0) {
            return VariantType.Reference;
        } else if (insertionLength == 1 && deletionLength == 1) {
            return VariantType.SNP;
        } else if (insertionLength == 0 && deletionLength > 0) {
            return VariantType.Deletion;
        } else if (deletionLength - insertionLength > 0 && new String(ref).endsWith(alts[ind - 1].getSeq().toString())) {
            //the second part will check if block substitutions possibly exist
            return VariantType.Deletion;
        } else if (insertionLength > 0 && deletionLength == 0) {
            return VariantType.Insertion;
        } else if (insertionLength == deletionLength) {
            return VariantType.MNP;
        }
        return VariantType.Complex;
    }


    /**
     * @return overall type of the variant considering all alleles
     */
    public VariantOverallType getType() {
        final Set<VariantType> variantTypes = EnumSet.noneOf(VariantType.class);
        if (compositions != null) {
            for (Variant c : compositions) {
                variantTypes.add(c.getType(c.getAllele(0)));
                variantTypes.add(c.getType(c.getAllele(1)));
            }
        } else {
            variantTypes.add(getType(getAllele(0)));
            variantTypes.add(getType(getAllele(1)));
        }

        if (variantTypes.size() == 1 && variantTypes.contains(VariantType.Reference)) {
            return VariantOverallType.Reference;
        }

        // The overall type depends on the kinds of non-reference alleles. If they are not the same, then return complex
        variantTypes.remove(VariantType.Reference);

        if (variantTypes.size() == 1) {
            if (variantTypes.contains(VariantType.SNP)) {
                return VariantOverallType.SNP;
            }
            if (variantTypes.contains(VariantType.Inversion)) {
                return VariantOverallType.Inversion;
            }
            if (variantTypes.contains(VariantType.Tandem_Duplication)) {
                return VariantOverallType.TandemDup;
            }
            if (variantTypes.contains(VariantType.Deletion)) {
                return VariantOverallType.Deletion;
            }
            if (variantTypes.contains(VariantType.Insertion)) {
                return VariantOverallType.Insertion;
            }
            if (variantTypes.contains(VariantType.Translocation_Deletion)) {
                return VariantOverallType.TransDel;
            }
            if (variantTypes.contains(VariantType.Translocation_Duplication)) {
                return VariantOverallType.TransDup;
            }
            if (variantTypes.contains(VariantType.Interspersed_Duplication)) {
                return VariantOverallType.InterDup;
            }
        } else if (variantTypes.contains(VariantType.Translocation_Deletion) && variantTypes.contains(VariantType.Translocation_Deletion)) {
            return VariantOverallType.Translocation;
        }

        /* Treat these as complex for now
        // check INDEL
        boolean is_indel = true;
        for (int a = 0; a < 2; a++) {
            if (allele[a] > 0 && !(getType(getAllele(a)) == Type.Deletion || getType(getAllele(a)) == Type.Insertion)) {
                is_indel = false;
                break;
            }
        }
        if (is_indel) {
            return VariantOverallType.INDEL;
        }

        // check MNP
        boolean is_mnp = true;
        for (int a = 0; a < 2; a++) {
            if (allele[a] > 0 && getType(getAllele(a)) != Type.MNP) {
                is_mnp = false;
                break;
            }
        }*/

        // otherwise it is complex since multiple kinds of variants occur at the same location
        return VariantOverallType.Complex;
    }

    /**
     * return alternative allele based on alternative allele index specified in GT field
     * alternative allele index = 1,2,...
     * @param ind
     * @return
     */
    public Alt getAlt(final int ind) {
        return (ind <= 0 || ind > alts.length) ? null : alts[ind - 1];
    }

    public void setAlt(final int ind, Alt alt) {
        if (ind > 0 && ind <= alts.length)
            alts[ind - 1] = alt;
    }

    public void setThreePrimeDistance(int threePrimeDistance) {
        this.threePrimeDistance = threePrimeDistance;
    }

    public void setFivePrimeDistance(int fivePrimeDistance) {
        this.fivePrimeDistance = fivePrimeDistance;
    }

    public void setLengthDifference(int lengthDifference) {
        this.lengthDifference = lengthDifference;
    }

    public String getFilter() {
        return filter;
    }

    public String getQual() {
        return qual;
    }

    public String getInfo() {
        return info;
    }

    public boolean isPhased() {
        return isPhased;
    }

    public boolean isRef() {
        return paternal == 0 && maternal == 0;
    }

    public String getVariantId() {
        return varId;
    }

    public byte[] getReference() {
        return ref;
    }

    public String getReferenceString() {
        try {
            return (refDeleted == null ? "" : refDeleted) + (ref == null ? "" : new String(ref, "US-ASCII"));
        } catch (UnsupportedEncodingException e) {
            e.printStackTrace();
            return "";
        }
    }

    public byte[] getRef_deleted() {
        return refDeleted.getBytes();
    }

    public int getCN(final int ind) {
        if (alts == null) {
          return 0;
        }
        return (ind <= 0 || ind > alts.length) ? 1 : alts[ind - 1].getCopyNumber();
    }

    /**
     * @return true if any of the alternate alleles has copy number greater than 1
     */
    public boolean hasCN() {
        if (alts == null) {
            return false;
        }
        boolean isCopyNumberPositive = false;
        for (Alt alt : alts) {
            if (alt.getCopyNumber() > 1) {
                isCopyNumberPositive = true;
            }
        }

        return isCopyNumberPositive;
    }

    /**
     *
     * @return . if ALT field is empty or unknown
     */
    public String alternativeAlleleString() {
        StringBuilder sbStr = new StringBuilder();
        String defaultResult = "N";
        if (alts == null) {
            return defaultResult;
        }
        for (int i = 0; i < alts.length; i++) {
            //if (i > 0 && alts[i].getSeq().toString().equals(alts[i - 1].toString())) {
                /*Marghoob suggested that two identical symbolic alternative alleles are
                allowed. so essentially we go back to original behavior of VarSim.
                 */
            if (i > 0) {
                sbStr.append(",");
            }
            if (alts[i].getSymbolicAllele() != null) {
                sbStr.append(alts[i].getSymbolicAllele().toString());
            } else if (alts[i].getSeq() != null) {
                if (alts[i].getSeq().isSeq()) {
                    sbStr.append(refDeleted).append(alts[i].getSeq().toString()).append(extraBase);
                } else {
                    sbStr.append(alts[i].getSeq().toString());
                }
            }
            sbStr.append(clippedSequence);
        }
        String result = sbStr.toString();
        if (result.length() == 0) {
            result = "<DEL>"; //most likely it's a deletion where all REF bases have been deleted (after trimming)
        }
        return result;
    }

    public int getNumberOfAlternativeAlleles() {
        return alts.length;
    }

    /**
     * Randomly swap the haploype
     */
    public void randomizeHaplotype() {
        if (rand == null) {
            log.error("Cannot randomize haplotype");
            log.error(toString());
            System.exit(1);
        }

        isPhased = true;
        if (rand.nextDouble() > 0.5) {
            return;
        }
        byte tmp = paternal;
        paternal = maternal;
        maternal = tmp;
    }

    /**
     * Randomize the genotype
     *
     * @param gender
     */
    public void randomizeGenotype(GenderType gender) {
        if (rand == null) {
            log.error("Cannot randomize genotype");
            log.error(toString());
            System.exit(1);
        }

        Genotypes g = new Genotypes(chr, gender, alts.length, rand);
        paternal = g.geno[0];
        maternal = g.geno[1];
        isPhased = true;
    }


    /**
    Tests if all of the alternate alleles with sequence are ACTGN
     */
    public boolean isAltACTGN() {
        for (Alt a : alts) {
            if (a.isSeq()) {
                if (!a.getSeq().toString().matches("[ACTGN]*")) {
                    return false;
                }
            }
        }
        return true;
    }

    /**
    Returns true if no genotype is missing
     and the variant is homozygous
     */
    public boolean isHom() {
        if (chr.isMT() && maternal < 0) {
            return false;
        }
        if (chr.isY() && paternal < 0) {
            return false;
        }
        if (maternal < 0 || paternal < 0) {
            return false;
        }
        return (paternal == maternal);
    }

    @Override
    public boolean equals(final Object o) {
        if (this == o) return true;
        if (!(o instanceof Variant)) return false;

        Variant variant = (Variant) o;

        if (referenceAlleleLength != variant.referenceAlleleLength) return false;
        if (isPhased != variant.isPhased) return false;
        if (maternal != variant.maternal) return false;
        if (paternal != variant.paternal) return false;
        if (pos != variant.pos) return false;
        if (wholeVariantIndex != variant.wholeVariantIndex) return false;
        if (splitVariantIndex != variant.splitVariantIndex) return false;
        if (!Arrays.equals(alts, variant.alts)) return false;
        if (chr != null ? !chr.equals(variant.chr) : variant.chr != null) return false;
        if (filter != null ? !filter.equals(variant.filter) : variant.filter != null) return false;
        if (info != null ? !info.equals(variant.info) : variant.info != null) return false;
        if (rand != null ? !rand.equals(variant.rand) : variant.rand != null) return false;
        if (!Arrays.equals(ref, variant.ref)) return false;
        if (refDeleted != null ? !refDeleted.equals(variant.refDeleted) : variant.refDeleted != null)
            return false;
        if (varId != null ? !varId.equals(variant.varId) : variant.varId != null) return false;
        if (originalType != variant.originalType) return false;

        return true;
    }

    @Override
    public int hashCode() {
        int result = hashCode;
        if (result != 0) return result;
       /*initial hash code cannot be zero, otherwise, we cannot
       tell the difference from hashCode if the first fields
       are zero.
        */
        result = rand != null ? rand.hashCode() : 17;
        result = 31 * result + splitVariantIndex;
        result = 31 * result + wholeVariantIndex;
        result = 31 * result + (originalType != null ? originalType.hashCode() : 0);
        result = 31 * result + pos;
        result = 31 * result + referenceAlleleLength;
        result = 31 * result + (ref != null ? Arrays.hashCode(ref) : 0);
        result = 31 * result + (chr != null ? chr.hashCode() : 0);
        result = 31 * result + (alts != null ? Arrays.hashCode(alts) : 0);
        result = 31 * result + (int) maternal;
        result = 31 * result + (int) paternal;
        result = 31 * result + (isPhased ? 1 : 0);
        result = 31 * result + (filter != null ? filter.hashCode() : 0);
        result = 31 * result + (info != null ? info.hashCode() : 0);
        result = 31 * result + (varId != null ? varId.hashCode() : 0);
        result = 31 * result + (refDeleted != null ? refDeleted.hashCode() : 0);
        hashCode = result;
        return result;
    }

    @Override
    public int compareTo(final Variant other) {
        final int chrCmp = chr.compareTo(other.chr);
        if (chrCmp != 0) {
            return chrCmp;
        }

        return getPos() - getRef_deleted().length - (other.getPos() - other.getRef_deleted().length);
    }

    public List<Integer> getSVLEN() {
        List<Integer> svlens = new ArrayList<Integer>();

        if (alts == null) {
            return svlens;
        }
        for (int i = 0; i < alts.length; i++) {
            VariantType t = getType(i + 1);
            int altLen = abs(alts[i].length());
            int svlen = 0;

            if (VariantType.Deletion.equals(t) || VariantType.Complex.equals(t) || VariantType.Insertion.equals(t)) {
                svlen = -referenceAlleleLength + altLen; // negative for deletions
            } else if (VariantType.MNP.equals(t) || VariantType.SNP.equals(t) || VariantType.Composite.equals(t)) {
                //composite variant contains more than 1 variant, no way to present as a regular variant
                svlen = 0;
            } else if (VariantType.Translocation_Deletion.equals(t)) {
                svlen = -abs(altLen);
            } else {
                //VariantType.Translocation_Duplication
                //VariantType.Inversion.equals(t)
//                VariantType.Interspersed_Duplication
//                VariantType.Tandem_Duplication
//                VariantType.Interspersed_Duplication
                svlen = altLen;
            }
            if (svlen != 0) {
                svlens.add(svlen);
            }
        }
        return svlens;
    }

    public String getLengthString() {
        List<Integer> svlens = getSVLEN();

        if (svlens.isEmpty()) {
            return "";
        } else {
            return "SVLEN=" + String.join(",", svlens.stream().map(s -> String.valueOf(s)).collect(Collectors.toList()));
        }
    }

    public void calculateExtraBase(final Sequence refSeq) {
        if (refDeleted.length() > 0) {
            return;
        }
        for (final Alt alt : alts) {
            if (alt.isSeq() && alt.length() == 0 && getPos() + referenceAlleleLength < refSeq.length()) {
                //why extrabase is only 1-bp long?
                extraBase = String.valueOf((char) refSeq.byteAt(getPos() + referenceAlleleLength));
                return;
            }
        }
        if (ref.length == 0) {
            extraBase = String.valueOf((char) refSeq.byteAt(getPos()));
        }
    }

    /**
     * @return a VCF record of the variant
     */
    @Override
    public String toString() {
        return toString(getGoodPaternal(), getGoodMaternal());
    }

    /**
     * @param paternal specified paternal allele
     * @param maternal specified maternal allele
     * @return the VCF record with prespecified genotype
     */
    public String toString(final int paternal, final int maternal) {
        StringBuilder sbStr = new StringBuilder();
        VariantOverallType t = getType();

        // chromosome name
        sbStr.append(chr == null ? "NA" : chr.toString());
        sbStr.append("\t");
        // start position
        sbStr.append(pos - (refDeleted == null ? 0 : refDeleted.length()));
        sbStr.append('\t');
        // variant id
        sbStr.append(varId);
        sbStr.append("\t");
        // ref allele
	String ref = getReferenceString() + extraBase;
      ref = ref + clippedSequence;
        if (ref.length() == 0) {
            ref = "N"; //make sure this field is not empty
        }
        sbStr.append(ref.toUpperCase());
        sbStr.append("\t");
        // alt alleles
        sbStr.append(alternativeAlleleString().toUpperCase());
        sbStr.append("\t");
        // variant quality
        sbStr.append(qual);
        sbStr.append("\t");
        // pass label
        sbStr.append(filter);
        sbStr.append("\t");
        // INFO
        List<String> infoFields = new ArrayList<String>();
        if (info != null) {
            // add original info fields before VarSim info fields
            // when duplicates are removed, original fields will be prioritized
            infoFields.addAll(Arrays.asList(info.split(";")));
        }
        infoFields.add("VARIANT_OVERALL_TYPE=" + getType());
        if (t == VariantOverallType.TandemDup) {
            infoFields.add("SVTYPE=DUP");
            infoFields.add(getLengthString());
            if(isInversed()) {
                infoFields.add("ISINV");
            }
        } else if (t == VariantOverallType.TransDup || t == VariantOverallType.InterDup) {
            infoFields.add("SVTYPE=DUP");
            if (getTraid() != null) {
                infoFields.add("TRAID=" + getTraid());
            }
            infoFields.add(getLengthString());
            //chr2,pos2,end2
            infoFields.add(String.format("CHR2=%s", StringUtilities.concatenateArray(getAllChr2(), ",")));
            infoFields.add(String.format("POS2=%s", StringUtilities.concatenateArray(getAllPos2(), ",")));
            infoFields.add(String.format("END2=%s", StringUtilities.concatenateArray(getAllEnd2(), ",")));
            if(isInversed()) {
                infoFields.add("ISINV");
            }
        } else if (t == VariantOverallType.Deletion) {
            infoFields.add("SVTYPE=DEL");
            infoFields.add(getLengthString());
        } else if (t == VariantOverallType.TransDel) {
            infoFields.add("SVTYPE=DEL");
            if (getTraid() != null) {
                infoFields.add("TRAID=" + getTraid());
            }
            infoFields.add(getLengthString());
        } else if (t == VariantOverallType.Inversion) {
            infoFields.add("SVTYPE=INV");
            infoFields.add(getLengthString());
        } else {
            infoFields.add(getLengthString());
        }
        /*
        checking -1 is for backward compatibility
        if output_distance_metric is disabled,
        then these INFO fields will not be set
        nor reported.
         */
        if (threePrimeDistance != -1) {
            infoFields.add(String.format("3PrimeDistance=%d", threePrimeDistance));
        }
        if (fivePrimeDistance != -1) {
            infoFields.add(String.format("5PrimeDistance=%d", fivePrimeDistance));
        }
        if (lengthDifference != -1) {
            infoFields.add(String.format("LengthDifference=%d", lengthDifference));
        }
	    sbStr.append(String.join(";", rmDupInfo(infoFields)));
        sbStr.append("\t");
        // label (GT)
        if (hasCN()) {
            sbStr.append("GT:CN\t");
        } else {
            sbStr.append("GT\t");
        }

        // for this one we need to work out which one is added
        if (paternal != -1 && maternal != -1) {
            sbStr.append(paternal);
            sbStr.append(this.isPhased() || paternal == maternal? "|" : "/");
            sbStr.append(maternal);
        } else if (paternal != -1){
            sbStr.append(paternal);
        } else if (maternal != -1) {
            sbStr.append(maternal);
        }

        if (hasCN()) {
            sbStr.append(":");
            if (paternal != -1 && maternal != -1) {
                sbStr.append(String.valueOf(getCN(paternal)));
                sbStr.append(this.isPhased() || paternal == maternal? "|" : "/");
                sbStr.append(String.valueOf(getCN(maternal)));
            } else if (paternal != -1) {
                sbStr.append(String.valueOf(getCN(paternal)));
            } else if (maternal != -1) {
                sbStr.append(String.valueOf(getCN(maternal)));
            }
        }
        return sbStr.toString().replaceAll(";+",";").replaceAll(";\t","\t").replaceAll("\t;","\t");
    }

    /**
     * remove duplicate and empty fields in info
     * first duplicate will be kept
     * duplicate is based on fieldname (before =), not on name=value pair
     * @param fields
     * @return
     */
    private List<String> rmDupInfo(List<String> fields) {
        Set<String> seen = new HashSet<>();
        List<String> uniqueNonEmptyFields = new ArrayList<>();
        for (String i : fields) {
            if (i.isEmpty() || i.equals(".")) {
                continue;
            } else {
                String fieldName = i.split("=")[0];
                if (!seen.contains(fieldName)) {
                    uniqueNonEmptyFields.add(i);
                    seen.add(fieldName);
                }
            }
        }
        return uniqueNonEmptyFields;
    }
    /**
     * write variant to a writer or write its compositions
     * to the writer using java8 streams
     *
     * @since 1.8
     * @param p
     */
    public void output(PrintWriter p) {
       if (compositions == null) {
           p.println(this);
       } else {
           compositions.stream().forEach(p::println);
       }
    }

    /**
     * this method returns a genotype number
     * -1 for unavailable genotype
     * 0 for reference allele
     * 1 for first alternative allele
     * 2,3,...
     *
     * when paternal genotype is
     * unknown, i.e. -1, then we try to return
     * something meaningful. however
     * this is by no means a very good return
     * value. proper fix involves explicitly handling
     * of missing genotype/phasing information.
     *
     * @return paternal genotype number
     */
    public byte getGoodPaternal() {
        if (paternal >= 0)
            return paternal; //not missing
        //missing
        //paternal genome does not have MT or X chromosome
        //works for male
        //TODO: for female, X is diploid, should not return -1
        if (chr.isMT() || chr.isX()) {
            return (byte) -1;
        }
        /*when Y or autosomal
        we try return a correct one
        however there is no guarantee that only one ALT is present
         */
        return (byte) (Math.min(1, alts.length));
    }

    /**
     * this method returns a genotype number
     * -1 for unavailable genotype
     * 0 for reference allele
     * 1 for first alternative allele
     * 2,3,...
     *
     * when maternal genotype is
     * unknown, i.e. -1, then we try to return
     * something meaningful. however
     * this is by no means a very good return
     * value. proper fix involves explicitly handling
     * of missing genotype/phasing information.
     *
     * @return maternal genotype number
     */
    public byte getGoodMaternal() {
        if (maternal >= 0)
            return maternal; //not missing
        //GT is missing, we try our best to return correct one
        //however there is no guarantee that only one ALT is specified
        //for a variant on MT
        if (chr.isMT()) {
          return (byte) (Math.min(1, alts.length));
        }
        //maternal genome does not have Y chromosome
        if (chr.isY()) {
            return (byte) -1;
        }
        /*chr=X or autosomal
        we try to return a correct one
        and make it work with getGoodPaternal(),i.e.
        getGoodPaternal() returns 1 by default,
        getGoodMaternal() returns 2 by default.
         */
        return (byte) (Math.min(2, alts.length));
    }

    public String getExtraBase() {
        return extraBase;
    }

    public Boolean isLengthImprecise() {
        return isLengthImprecise;
    }

    public void setLengthImprecise(Boolean lengthImprecise) {
        isLengthImprecise = lengthImprecise;
    }
}
