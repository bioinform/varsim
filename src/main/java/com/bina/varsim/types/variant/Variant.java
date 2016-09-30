package com.bina.varsim.types.variant;

//--- Java imports ---

import com.bina.intervalTree.SimpleInterval1D;
import com.bina.varsim.types.*;
import com.bina.varsim.util.SimpleReference;
import org.apache.log4j.Logger;

import java.io.UnsupportedEncodingException;
import java.util.Arrays;
import java.util.Random;

public class Variant implements Comparable<Variant>{
    private final static Logger log = Logger.getLogger(Variant.class.getName());
    public int idx = 0; // this is hopefully a unique index, for the split variants
    public int full_idx = 0; // this is hopefully a unique index, for the whole variants
    // this is the type before the variant was split into canonical ones
    public VariantOverallType original_type = null;
    // use a seed for reproducibility, should be an option or global
    private Random _rand = null;
    //what does _del mean exactly? length of reference allele?
    private int _pos = -1, _del = -1;
    private byte[] _ref;
    private ChrString _chr;
    private FlexSeq[] _alts;
    private byte _maternal = 0, _paternal = 0; // -1 for not avaliable
    private boolean _isPhased = false; // Phasing
    private String _filter;
    private String _var_id;
    // this is when the reference base is deleted
    // if it is the same as the first alt base
    //so _ref_deleted.length() <= 1 is always true?
    private String _ref_deleted;
    private String extraBase = "";

    public Variant(Random rand) {
        // TODO define some methods to determine if a Variant is uninitialised
        _rand = rand;
    }

    public Variant(ChrString chr, int pos, int del, byte[] ref,
                   FlexSeq[] alts, byte[] phase, boolean isPhased, String var_id, String filter,
                   String ref_deleted) {
        this(chr, pos, del, ref, alts, phase, isPhased, var_id, filter, ref_deleted, null);
    }

    public Variant(ChrString chr, int pos, int del, byte[] ref,
                   FlexSeq[] alts, byte[] phase, boolean isPhased, String var_id, String filter,
                   String ref_deleted, Random rand) {
        this(rand);

        _filter = filter;
        _var_id = var_id;

        _chr = chr;
        _pos = pos;
        _del = del;

        // TODO we should put the reference matching code here
        _ref = ref.clone();
        _ref_deleted = ref_deleted;
        _alts = new FlexSeq[alts.length];
        for (int i = 0; i < alts.length; i++) {
            if (alts[i] != null) {
                _alts[i] = new FlexSeq(alts[i]);
            } else {
                _alts[i] = null;
            }
        }

        _paternal = phase[0];
        _maternal = phase[1];
        _isPhased = isPhased;
    }

    public Variant(final Variant var) {
        _filter = var._filter;
        _var_id = var._var_id;
        _chr = var._chr;
        _pos = var._pos;
        _del = var._del;
        _ref = var._ref.clone();
        _ref_deleted = var._ref_deleted;
        _alts = new FlexSeq[var._alts.length];
        for (int i = 0; i < var._alts.length; i++) {
            if (var._alts[i] != null) {
                _alts[i] = new FlexSeq(var._alts[i]);
            } else {
                _alts[i] = null;
            }
        }

        _paternal = var._paternal;
        _maternal = var._maternal;
        _isPhased = var._isPhased;
        _rand = var._rand;
    }

    /**
     * @return Chromosome variant is on
     */
    public ChrString getChr() {
        return _chr;
    }

    /**
     * @return Start position of variant
     */
    public int getPos() {
        return _pos;
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
    public boolean setNovelPosition(int pos, SimpleReference ref) {

        // replace ref
        int len = _ref.length;

        if (len > 0) {
            byte[] temp_ref = ref.byteRange(_chr, pos, pos + len);

            for (byte b : temp_ref) {
                if (b == 'N') {
                    // don't allow N's
                    log.warn("N found at " + pos + " to " + (pos + len));
                    return false;
                }
            }

            for (FlexSeq f : _alts) {
                if (f.getSeq() != null) {
                    // make sure there is no prefix the same
                    for (int i = 0; i < temp_ref.length; i++) {
                        if (i < f.getSeq().length) {
                            if (temp_ref[i] == f.getSeq()[i]) {
                                log.warn("Same ref at alt at " + pos + " to " + (pos + len));
                                return false;
                            } else {
                                break;
                            }
                        }
                    }

                }
            }

            _ref = temp_ref;
        }

        // replace ref_deleted
        len = _ref_deleted.length();
        if (len > 0) {
            try {
                byte[] deleted_temp = ref.byteRange(_chr, pos - len, pos);
                _ref_deleted = new String(deleted_temp, "US-ASCII");
            } catch (UnsupportedEncodingException e) {
                e.printStackTrace();
            }
        }

        // do the replacing
        _pos = pos;

        return true;
    }

    /**
     * @param id variant id. usually the dbSNP id
     */
    public void setVarID(String id) {
        _var_id = id;
    }

    /**
     * return the length of the deletion.
     * This is currently simply the length of the reference sequence that will
     * be replaced
     *
     * @return the length of the deletion
     */
    public int deletion() {
        return _del;
    }

    /**
     * @return maternal allele index, 0 if reference
     */
    public int maternal() {
        return _maternal;
    }

    /**
     * @return paternal allele index, 0 if reference
     */
    public int paternal() {
        return _paternal;
    }

    /**
     * @param ind index of allele
     * @return the insertion sequence as a string
     */
    public byte[] insertion(int ind) {
        if (ind <= 0 || ind > _alts.length)
            return null;
        return _alts[ind - 1].getSeq();
    }

    /**
     * The length of an alternate allele, this is usually an insertion sequence
     * But in the case of SVs, not necessarily
     *
     * @param ind index of allele
     * @return the length of that allele
     */
    public int insertion_len(int ind) {
        if (ind <= 0 || ind > _alts.length)
            return 0;
        return _alts[ind - 1].length();
    }

    // if it is a simple indel, it is just the length
    // if it is a complex variant, this is the maximum length of the insertion
    // and deletion
    public int maxLen(int ind) {
        if (ind <= 0 || ind > _alts.length)
            return 0;
        return Math.max(_del, _alts[ind - 1].length());
    }

    public int maxLen() {
        int len = 0;
        for (int i = 0; i < 2; i++) {
            len = Math.max(maxLen(get_allele(i)), len);
        }
        return len;
    }

    // this is the minimum length of the variants
    public int minLen() {
        int len = maxLen(get_allele(0));
        len = Math.min(maxLen(get_allele(1)), len);

        return len;
    }

    /*
    gets the interval enclosing the variant on the reference genome
    */
    public SimpleInterval1D get_interval(int ind) {
        if (ind == 0 || _del == 0) {
            return new SimpleInterval1D(_pos, _pos);
        }

        return new SimpleInterval1D(_pos, _pos + _del - 1);
    }

    /*
    gets the interval for the variant, accounting for variant size
    */
    public SimpleInterval1D get_var_interval(int ind) {
        try {
            if (ind == 0) {
                return new SimpleInterval1D(_pos, _pos);
            }

            // TODO hmm unsafe
            if (maxLen(ind) == Integer.MAX_VALUE) {
                return new SimpleInterval1D(_pos, _pos);
            } else {
                return new SimpleInterval1D(_pos, _pos + maxLen(ind) - 1);
            }
        } catch (RuntimeException e) {
            log.error("Bad variant interval: " + toString());
            log.error("_pos: " + _pos);
            log.error("ind: " + ind);
            log.error("maxLen(ind): " + maxLen(ind));
            e.printStackTrace();
            return null;
        }
    }

    // union of intervals from the genotypes
    public SimpleInterval1D get_geno_interval() {
        return get_interval(getgood_paternal()).union(get_interval(getgood_maternal()));
    }

    public SimpleInterval1D get_geno_var_interval() {
        return get_var_interval(getgood_paternal()).union(get_var_interval(getgood_maternal()));
    }

    public Genotypes getGeno() {
        return new Genotypes(getgood_paternal(), getgood_maternal());
    }

    /*
    * 0 = paternal
    * 1 = maternal
    * otherwise returns -1
     */
    public int get_allele(int parent) {
        if (parent == 0) {
            return getgood_paternal();
        } else if (parent == 1) {
            return getgood_maternal();
        }
        return -1;
    }

    public void set_allele(int parent, byte allele) {
        if (parent == 0) {
            _paternal = allele;
        } else if (parent == 1) {
            _maternal = allele;
        }
    }

    // TODO this is wrong, but it only effects the count of variant bases
    public int variantBases() {
        int ret = _del;
        for (FlexSeq _alt : _alts) {
            if (_del != _alt.length()) {
                ret += _alt.length();
            }
        }
        return ret;
    }

    /**
     * @param ind index of allele (starts at 1, 0 is reference)
     * @return type of allele at index ind
     */
    public VariantType getType(int ind) {
        if (ind == 0) {
            return VariantType.Reference;
        }

        /*when variant type is explicitly specified in VCF,
        alt allele will be set to the correct type,
        we can be assured that correct variant type is
        returned.
         */
        FlexSeq.Type type = _alts[ind - 1].getType();
        switch (type) {
            case DUP:
                return VariantType.Tandem_Duplication;
            case INS:
                return VariantType.Insertion;
            case INV:
                return VariantType.Inversion;
            default:
                break;
        }

        int inslen = insertion_len(ind);
        int dellen = _del;
        if (inslen == 0 && dellen == 0) {
            return VariantType.Reference;
        } else if (inslen == 1 && dellen == 1) {
            return VariantType.SNP;
        } else if (inslen == 0 && dellen > 0) {
            return VariantType.Deletion;
        } else if (dellen - inslen > 0 && new String(_ref).endsWith(_alts[ind - 1].toString())) {
            return VariantType.Deletion;
        } else if (inslen > 0 && dellen == 0) {
            return VariantType.Insertion;
        } else if (inslen == dellen) {
            return VariantType.MNP;
        }
        return VariantType.Complex;
    }

    /**
     * @return overall type of the variant considering both alleles
     */
    public VariantOverallType getType() {
        int[] allele = {get_allele(0), get_allele(1)};

        // check Reference
        boolean is_ref = true;
        for (int a = 0; a < 2; a++) {
            if (getType(allele[a]) != VariantType.Reference) {
                is_ref = false;
                break;
            }
        }
        if (is_ref) {
            return VariantOverallType.Reference;
        }

        // check SNP
        boolean is_snp = true;
        for (int a = 0; a < 2; a++) {
            if (allele[a] > 0 && getType(allele[a]) != VariantType.SNP) {
                is_snp = false;
                break;
            }
        }
        if (is_snp) {
            return VariantOverallType.SNP;
        }

        // check INV
        boolean is_inv = true;
        for (int a = 0; a < 2; a++) {
            if (allele[a] > 0 && getType(allele[a]) != VariantType.Inversion) {
                is_inv = false;
                break;
            }
        }
        if (is_inv) {
            return VariantOverallType.Inversion;
        }

        // check DUP
        boolean is_dup = true;
        for (int a = 0; a < 2; a++) {
            if (allele[a] > 0 && getType(allele[a]) != VariantType.Tandem_Duplication) {
                is_dup = false;
                break;
            }
        }
        if (is_dup) {
            return VariantOverallType.Tandem_Duplication;
        }

        // check Deletion
        boolean is_del = true;
        for (int a = 0; a < 2; a++) {
            if (allele[a] > 0 && getType(allele[a]) != VariantType.Deletion) {
                is_del = false;
                break;
            }
        }
        if (is_del) {
            return VariantOverallType.Deletion;
        }

        // check DUP
        boolean is_ins = true;
        for (int a = 0; a < 2; a++) {
            if (allele[a] > 0 && getType(allele[a]) != VariantType.Insertion) {
                is_ins = false;
                break;
            }
        }
        if (is_ins) {
            return VariantOverallType.Insertion;
        }

        /* Treat these as complex for now
        // check INDEL
        boolean is_indel = true;
        for (int a = 0; a < 2; a++) {
            if (allele[a] > 0 && !(getType(get_allele(a)) == Type.Deletion || getType(get_allele(a)) == Type.Insertion)) {
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
            if (allele[a] > 0 && getType(get_allele(a)) != Type.MNP) {
                is_mnp = false;
                break;
            }
        }
        if (is_mnp) {
            return VariantOverallType.MNP;
        }
        */

        // otherwise it is complex
        return VariantOverallType.Complex;
    }

    public FlexSeq getAlt(int ind) {
        if (ind <= 0 || ind > _alts.length)
            return null;
        return _alts[ind - 1];
    }

    public void setAlt(int ind, FlexSeq alt) {
        if (ind <= 0 || ind > _alts.length) {
            return;
        }
        _alts[ind - 1] = alt;
    }

    public String getFilter() {
        return _filter;
    }

    public boolean isPhased() {
        return _isPhased;
    }

    public boolean isRef() {
        return _paternal == 0 && _maternal == 0;
    }

    public String getVar_id() {
        return _var_id;
    }

    public byte[] getRef() {
        return _ref;
    }

    public String getOrig_Ref() {
        try {
            return _ref_deleted + new String(_ref, "US-ASCII");
        } catch (UnsupportedEncodingException e) {
            e.printStackTrace();
            return "";
        }
    }

    public String getRef_deleted() {
        return _ref_deleted;
    }

    public int getCN(int ind) {
        if (ind <= 0 || ind > _alts.length)
            return 0;
        return _alts[ind - 1].getCopy_num();
    }

    /**
     * @return true if any of the alternate alleles has copy number greater than 1
     */
    public boolean hasCN() {
        boolean CN_positive = false;
        for (FlexSeq _alt : _alts) {
            if (_alt.getCopy_num() > 1) {
                CN_positive = true;
            }
        }

        return CN_positive;
    }

    public StringBuilder alt_string() {
        StringBuilder sbStr = new StringBuilder();
        for (int i = 0; i < _alts.length; i++) {
            if (i > 0) {
                sbStr.append(",");
            }
            if (_alts[i].isSeq()) {
                sbStr.append(_ref_deleted).append(_alts[i].toString()).append(extraBase);
            } else {
                sbStr.append(_alts[i].toString());
            }
        }
        return sbStr;
    }

    public int get_num_alt() {
        return _alts.length;
    }

    /**
     * Randomly swap the haploype
     */
    public void randomizeHaplotype() {
        if (_rand == null) {
            log.error("Cannot randomize haplotype");
            log.error(toString());
            System.exit(1);
        }

        if (_rand.nextDouble() > 0.5) {
            return;
        }
        byte tmp = _paternal;
        _paternal = _maternal;
        _maternal = tmp;
    }

    /**
     * Randomize the genotype
     *
     * @param gender
     */
    public void randomizeGenotype(GenderType gender) {
        if (_rand == null) {
            log.error("Cannot randomize genotype");
            log.error(toString());
            System.exit(1);
        }

        Genotypes g = new Genotypes(_chr, gender, _alts.length, _rand);
        _paternal = g.geno[0];
        _maternal = g.geno[1];
    }


    /*
    Tests if all of the alternate alleles with sequence are ACTGN
     */
    public boolean isAltACTGN() {
        for (FlexSeq a : _alts) {
            if (a.isSeq()) {
                if (!a.toString().matches("[ACTGN]*")) {
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
        return (_paternal == _maternal);
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;

        Variant variant = (Variant) o;

        if (_del != variant._del) return false;
        if (_isPhased != variant._isPhased) return false;
        if (_maternal != variant._maternal) return false;
        if (_paternal != variant._paternal) return false;
        if (_pos != variant._pos) return false;
        if (full_idx != variant.full_idx) return false;
        if (idx != variant.idx) return false;
        if (!Arrays.equals(_alts, variant._alts)) return false;
        if (_chr != null ? !_chr.equals(variant._chr) : variant._chr != null) return false;
        if (_filter != null ? !_filter.equals(variant._filter) : variant._filter != null) return false;
        if (_rand != null ? !_rand.equals(variant._rand) : variant._rand != null) return false;
        if (!Arrays.equals(_ref, variant._ref)) return false;
        if (_ref_deleted != null ? !_ref_deleted.equals(variant._ref_deleted) : variant._ref_deleted != null)
            return false;
        if (_var_id != null ? !_var_id.equals(variant._var_id) : variant._var_id != null) return false;
        if (original_type != variant.original_type) return false;

        return true;
    }

    @Override
    public int hashCode() {
        int result = _rand != null ? _rand.hashCode() : 0;
        result = 31 * result + idx;
        result = 31 * result + full_idx;
        result = 31 * result + (original_type != null ? original_type.hashCode() : 0);
        result = 31 * result + _pos;
        result = 31 * result + _del;
        result = 31 * result + (_ref != null ? Arrays.hashCode(_ref) : 0);
        result = 31 * result + (_chr != null ? _chr.hashCode() : 0);
        result = 31 * result + (_alts != null ? Arrays.hashCode(_alts) : 0);
        result = 31 * result + (int) _maternal;
        result = 31 * result + (int) _paternal;
        result = 31 * result + (_isPhased ? 1 : 0);
        result = 31 * result + (_filter != null ? _filter.hashCode() : 0);
        result = 31 * result + (_var_id != null ? _var_id.hashCode() : 0);
        result = 31 * result + (_ref_deleted != null ? _ref_deleted.hashCode() : 0);
        return result;
    }

    @Override
    public int compareTo(Variant other) {
        final int chrCmp = _chr.compareTo(other._chr);
        if (chrCmp != 0) {
            return chrCmp;
        }

        return getPos() - getRef_deleted().length() - (other.getPos() - other.getRef_deleted().length());
    }

    public String getLength() {
        StringBuilder len = new StringBuilder();

        for (int i = 0; i < _alts.length; i++) {
            if (i > 0) {
                len.append(',');
            }

            VariantType t = getType(i + 1);

            if (t == VariantType.Deletion) {
                len.append(-_del + _alts[i].length()); // negative for deletions
            } else if (t == VariantType.Complex) {
                int alt_len = _alts[i].length();
                if (_del > alt_len) {
                    len.append(-_del);
                } else {
                    len.append(alt_len);
                }
            } else {
                len.append(_alts[i].length());
            }
        }

        return len.toString();
    }

    public void calculateExtraBase(final Sequence refSeq) {
        for (final FlexSeq alt : _alts) {
            if (alt.isSeq() && alt.length() == 0 && getPos() + _del < refSeq.length()) {
                //why extrabase is only 1-bp long?
                extraBase = String.valueOf((char) refSeq.byteAt(getPos() + _del ));
            }
        }
    }


    /**
     * @param sbStr will build a VCF record without genotype
     */
    // TODO, this should be self contained and output a VCF record
    private void buildVCFstr(StringBuilder sbStr) {
        // chromosome name
        sbStr.append(_chr.toString());
        sbStr.append("\t");
        // start position
        sbStr.append(_pos - _ref_deleted.length());
        sbStr.append('\t');
        // variant id
        sbStr.append(_var_id);
        sbStr.append("\t");
        // ref allele
        sbStr.append(getOrig_Ref());
        sbStr.append("\t");
        // alt alleles
        sbStr.append(alt_string().toString());
        sbStr.append("\t");
        // variant quality
        sbStr.append(".\t");
        // pass label
        sbStr.append(_filter);
        sbStr.append("\t");
        // INFO
        if (getType() == VariantOverallType.Tandem_Duplication) {
            sbStr.append("SVTYPE=DUP;");
            sbStr.append("SVLEN=");
            sbStr.append(getLength());
        } else if (getType() == VariantOverallType.Inversion) {
            sbStr.append("SVTYPE=INV;");
            sbStr.append("SVLEN=");
            sbStr.append(getLength());
        } else {
            sbStr.append("SVLEN=");
            sbStr.append(getLength());
        }
        sbStr.append("\t");

        // label (GT)
        if (hasCN()) {
            sbStr.append("CN:GT\t");
        } else {
            sbStr.append("GT\t");
        }

        if (hasCN()) {
            sbStr.append(String.valueOf(getCN(getgood_paternal())));
            sbStr.append("|");
            sbStr.append(String.valueOf(getCN(getgood_maternal())));
            sbStr.append(":");
        }

    }

    /**
     * @return a VCF record of the variant
     */
    public String toString() {
        StringBuilder sbStr = new StringBuilder();

        buildVCFstr(sbStr);


        sbStr.append(getgood_paternal());
        sbStr.append("|");
        sbStr.append(getgood_maternal());

        return sbStr.toString();
    }

    /**
     * @param paternal specified paternal allele
     * @param maternal specified maternal allele
     * @return the VCF record with prespecified genotype
     */
    public String toString(int paternal, int maternal) {
        StringBuilder sbStr = new StringBuilder();

        buildVCFstr(sbStr);

        // for this one we need to work out which one is added
        sbStr.append(paternal);
        sbStr.append("|");
        sbStr.append(maternal);

        return sbStr.toString();
    }

    public byte getgood_paternal() {
        if (_paternal < 0) {
            return 1;
        } else {
            return _paternal;
        }
    }

    public byte getgood_maternal() {
        if (_maternal < 0) {
            return 1;
        } else {
            return _maternal;
        }
    }

    public String getExtraBase() {
        return extraBase;
    }

}
