package org.go;

public class Gene
{

    private String symbol;
    private double fc;
    private boolean sig;

    public Gene(String id, double fc, boolean sig)
    {
        this.symbol = id;
        this.fc = fc;
        this.sig = sig;
    }

    public String getSymbol()
    {
        return symbol;
    }

    public void setSymbol(String symbol)
    {
        this.symbol = symbol;
    }

    public double getFc()
    {
        return fc;
    }

    public void setFc(double fc)
    {
        this.fc = fc;
    }

    public boolean isSig()
    {
        return sig;
    }

    public void setSig(boolean sig)
    {
        this.sig = sig;
    }
}
