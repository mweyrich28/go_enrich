package org.go;

import java.util.ArrayList;

public class GOEntry
{

    private String id;
    private String annot;
    private ArrayList<GOEntry> parents;
    private ArrayList<GOEntry> children;
    private ArrayList<Gene> genes;

    public GOEntry(String id, String annot)
    {
        this.id = id;
        this.annot = annot;
        this.genes = new ArrayList<>();
//        this.parents = new ArrayList<>();
        this.children = new ArrayList<>();
    }

    public void setAnnot(String annot)
    {
        if (this.annot == null)
        {
            this.annot = annot;
        }
    }

    public String getId()
    {
        return id;
    }

    public void addGene(Gene gene)
    {
        this.genes.add(gene);
    }

    public void addParents(ArrayList<String> parents)
    {

    }

    public void addChild(GOEntry child)
    {
        this.children.add(child);
    }

    public ArrayList<GOEntry> getChildren()
    {
        return children;
    }

    @Override
    public String toString()
    {
        return this.id + "\t" + this.annot + "\n\t" + children.size();
    }

}
