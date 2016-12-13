package com.bina.varsim.types.variant;

//--- Java imports ---

import com.bina.intervalTree.SimpleInterval1D;
import com.bina.varsim.types.*;
import com.bina.varsim.types.variant.alt.Alt;
import com.bina.varsim.util.SimpleReference;
import org.apache.log4j.Logger;

import java.io.UnsupportedEncodingException;
import java.util.*;

import static java.lang.Character.getType;

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
    private String varId;
    // this is when the reference base is deleted
    // if it is the same as the first alt base
    //so refDeleted.length() <= 1 is always true?
    private String refDeleted;
    private String extraBase = "";
    private ChrString[] chr2;
    private int[] pos2;
    private int[] end2;
    private int end;
    private Boolean isinv; //is sequence inverted? useful for interspersed dup, translocation dup
    private String traid; //translocation ID
    private volatile int hashCode; //caching hashcode
    private List<Variant> compositions; //when this variant is a composite variant, this field will store constitutional variants


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
     * @param idx
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
        private String varId;
        // this n the reference base is deleted
        // if ite same as the first alt base
        //so refd.length() <= 1 is always true?
        private String refDeleted;
        private String extraBase = "";
        private ChrString[] chr2;
        private int[] pos2;
        private int[] end2;
        private int end;
        private Boolean isinv; //is sequence inverted? useful for interspersed dup, translocation dup
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
        public Builder end(final int end) {
            this.end = end;
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
        public Builder traid(final String id) {
            this.traid = id;
            return this;
        }
        public Builder compositions(final List<Variant> compositions) {
            this.compositions = compositions;
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
        this.chr = builder.chr;
        this.pos = builder.pos;
        this.referenceAlleleLength = builder.referenceAlleleLength;

        this.ref = builder.ref;
        this.refDeleted = builder.refDeleted;
        this.alts = builder.alts;
        this.chr2 = builder.chr2;
        this.pos2 = builder.pos2;
        this.end2 = builder.end2;
        this.paternal = builder.paternal;
        this.maternal = builder.maternal;
        this.isPhased = builder.isPhased;
        this.end = builder.end;
        this.traid = builder.traid;
        this.isinv = builder.isinv;
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
        varId = var.varId;
        chr = var.chr;
        pos = var.pos;
        referenceAlleleLength = var.referenceAlleleLength;
        ref = var.ref == null ? null : var.ref.clone();
        refDeleted = var.refDeleted;
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
        this.end = var.end;
        paternal = var.paternal;
        maternal = var.maternal;
        isPhased = var.isPhased;
        rand = var.rand;
        isinv = var.isinv;
        traid = var.traid;
        compositions = var.getCompositions();
    }

    /**
     * @return Chromosome variant is on
     */
    public ChrString getChr() {
        if (chr == null) {
            throw new UnsupportedOperationException("ERROR: no legitimate chromosome name available!");
        }
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
        return this.chr2;
    }

    public int[] getAllPos2() {
        return this.pos2;
    }

    public int[] getAllEnd2() {
        return this.end2;
    }

    public int getEnd() {
        return this.end;
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

    // union of intervals from the genotypes
    public SimpleInterval1D getGenotypeUnionAlternativeInterval() {
        return getAlternativeAlleleInterval(getGoodPaternal()).union(getAlternativeAlleleInterval(getGoodMaternal()));
    }

    public SimpleInterval1D getGenotypeUnionVariantInterval() {
        return getVariantInterval(getGoodPaternal()).union(getVariantInterval(getGoodMaternal()));
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
        if (ind == 0) {
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
        Alt alt = alts[ind - 1];
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
            return VariantOverallType.CutPasteTrans;
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

    public String getFilter() {
        return filter;
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
            return refDeleted + new String(ref, "US-ASCII");
        } catch (UnsupportedEncodingException e) {
            e.printStackTrace();
            return "";
        }
    }

    public byte[] getRef_deleted() {
        return refDeleted.getBytes();
    }

    public int getCN(final int ind) {
        return (ind <= 0 || ind > alts.length) ? 0 : alts[ind - 1].getCopyNumber();
    }

    /**
     * @return true if any of the alternate alleles has copy number greater than 1
     */
    public boolean hasCN() {
        boolean isCopyNumberPositive = false;
        for (Alt alt : alts) {
            if (alt.getCopyNumber() > 1) {
                isCopyNumberPositive = true;
            }
        }

        return isCopyNumberPositive;
    }

    public String alternativeAlleleString() {
        StringBuilder sbStr = new StringBuilder();
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
        }
        return sbStr.toString();
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
    }


    /*
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

    /*
    Returns true if the variant is homozygous
     */
    public boolean isHom() {
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

    public String getLengthString() {
        StringBuilder len = new StringBuilder();

        for (int i = 0; i < alts.length; i++) {
            if (i > 0) {
                len.append(',');
            }

            VariantType t = getType(i + 1);

            if (VariantType.Deletion.equals(t)) {
                len.append(-referenceAlleleLength + alts[i].length()); // negative for deletions
            } else if (VariantType.Translocation_Deletion.equals(t)) {
                len.append(-referenceAlleleLength); // negative for deletions
            } else if (VariantType.Complex.equals(t)) {
                int alt_len = alts[i].length();
                if (referenceAlleleLength > alt_len) {
                    len.append(-referenceAlleleLength);
                } else {
                    len.append(alt_len);
                }
            } else {
                len.append(alts[i].length());
            }
        }

        return len.toString();
    }

    public void calculateExtraBase(final Sequence refSeq) {
        for (final Alt alt : alts) {
            if (alt.isSeq() && alt.length() == 0 && getPos() + referenceAlleleLength < refSeq.length()) {
                //why extrabase is only 1-bp long?
                extraBase = String.valueOf((char) refSeq.byteAt(getPos() + referenceAlleleLength));
            }
        }
    }


    /**
     * @param sbStr will build a VCF record without genotype
     */
    // TODO, this should be self contained and output a VCF record
    private void buildVCFstr(final StringBuilder sbStr) {
        // chromosome name
        sbStr.append(chr.toString());
        sbStr.append("\t");
        // start position
        sbStr.append(pos - refDeleted.length());
        sbStr.append('\t');
        // variant id
        sbStr.append(varId);
        sbStr.append("\t");
        // ref allele
        sbStr.append(getReferenceString());
        sbStr.append("\t");
        // alt alleles
        sbStr.append(alternativeAlleleString());
        sbStr.append("\t");
        // variant quality
        sbStr.append(".\t");
        // pass label
        sbStr.append(filter);
        sbStr.append("\t");
        // INFO
        if (getType() == VariantOverallType.TandemDup || getType() == VariantOverallType.InterDup || getType() == VariantOverallType.TransDup) {
            sbStr.append("SVTYPE=DUP;");
            if (traid != null)
                sbStr.append("TRAID=" + traid + ";");
            sbStr.append("SVLEN=");
            sbStr.append(getLengthString());
        } else if (getType() == VariantOverallType.Deletion || getType() == VariantOverallType.TransDel) {
            sbStr.append("SVTYPE=DEL;");
            if (traid != null)
                sbStr.append("TRAID=" + traid + ";");
            sbStr.append("SVLEN=");
            sbStr.append(getLengthString());
        } else if (getType() == VariantOverallType.Inversion) {
            sbStr.append("SVTYPE=INV;");
            sbStr.append("SVLEN=");
            sbStr.append(getLengthString());
        } else {
            sbStr.append("SVLEN=");
            sbStr.append(getLengthString());
        }
        sbStr.append("\t");

        // label (GT)
        if (hasCN()) {
            sbStr.append("CN:GT\t");
        } else {
            sbStr.append("GT\t");
        }

        if (hasCN()) {
            sbStr.append(String.valueOf(getCN(getGoodPaternal())));
            sbStr.append("|");
            sbStr.append(String.valueOf(getCN(getGoodMaternal())));
            sbStr.append(":");
        }

    }

    /**
     * @return a VCF record of the variant
     */
    public String toString() {
        StringBuilder sbStr = new StringBuilder();

        buildVCFstr(sbStr);


        sbStr.append(getGoodPaternal());
        sbStr.append("|");
        sbStr.append(getGoodMaternal());

        return sbStr.toString();
    }

    /**
     * @param paternal specified paternal allele
     * @param maternal specified maternal allele
     * @return the VCF record with prespecified genotype
     */
    public String toString(final int paternal, final int maternal) {
        StringBuilder sbStr = new StringBuilder();

        buildVCFstr(sbStr);

        // for this one we need to work out which one is added
        sbStr.append(paternal);
        sbStr.append("|");
        sbStr.append(maternal);

        return sbStr.toString();
    }

    /**
     * this method returns a genotype
     * however, why does it return 1
     * when paternal genotype is < 0?
     * 0 for reference allele
     * 1 for first alternative allele
     * -1 for unavailable genotype
     *
     * here it guarantees to return
     * genotype 1 or 2
     *
     * @return
     */
    public byte getGoodPaternal() {
        return (paternal < 0) ? 1 : paternal;

    }

    public byte getGoodMaternal() {
        return (maternal < 0) ? 1: maternal;
    }

    public String getExtraBase() {
        return extraBase;
    }

}
