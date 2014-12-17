package com.binatechnologies.varsim;

import org.apache.commons.lang3.StringUtils;
import org.apache.log4j.Logger;

import java.util.ArrayList;

/**
 *
 * This current implementation allows for duplicates
 *
 * Created by johnmu on 12/6/14.
 */
public class IntervalTree<Key extends Interval1D> {
    private final static Logger log = Logger.getLogger(IntervalTree.class.getName());
    IntervalTreeNode<Key> root = null;

    IntervalTree(){}

    /**
     * Adds key to the root
     * @param k
     */
    void add(Key k){
        add(root,k);
    }

    /**
     * Adds key to a particular node in the tree
     * @param head
     * @param k
     */
    void add(IntervalTreeNode<Key> head, Key k){
        if(head == null){
            // Uninitialised node, simply add the key
            head = new IntervalTreeNode<Key>(k.getCenter());
            head.addKey(k);
            return;
        }else{
            // Try to add the key
            int compVal = head.addKey(k);
            if(compVal == 0){
                // added to center, we can stop now :D
                return;
            }else if(compVal < 0){
                add(head.getLeft(), k);
            }else{
                add(head.getRight(), k);
            }
        }
    }


    /**
     * Gets keys overlapping k from the root
     * @param k
     * @return
     */
    ArrayList<Key> getContains(Key k){
        return getContains(root,k);
    }

    /**
     *
     * There are some optimizations here if the left and right trees are augmented with more information
     *
     * @param k key to search for
     * @return List of Keys that overlap the provided key
     */
    ArrayList<Key> getContains(IntervalTreeNode<Key> head, Key k){
        ArrayList<Key> retVal = head.getOverlaps(k);
        int compVal = head.checkKey(k);
        if(compVal == 0){
            retVal.addAll(getContains(head.getLeft(),k));
            retVal.addAll(getContains(head.getRight(),k));
        }else if(compVal < 0){
            retVal.addAll(getContains(head.getLeft(),k));
        }else{
            retVal.addAll(getContains(head.getRight(),k));
        }
        return retVal;
    }

    boolean contains(){
        // stub
        return false;
    }

    void rotateRight(){
        // stub
    }

    void rotateLeft(){
        // stub
    }

    public static void main(String args[]){
        IntervalTree<Interval1D> test_tree = new IntervalTree<Interval1D>();
        test_tree.add(new Interval1D(10,20));
        ArrayList<Interval1D> output = test_tree.getContains(new Interval1D(10,10));
        System.out.println(StringUtils.join(output, ','));
    }

}
