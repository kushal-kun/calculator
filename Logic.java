
import java.util.*;
public class Logic {
    static Scanner sc =  new Scanner(System.in);
    public static void main(String[] args) {
       
    }
    public static double nMean(double values[]){
        double sum = 0;
        for(int i=0;i<values.length;i++){
            sum+=values[i];
        }
        return sum/values.length;
    }
    public static double gMean(double values[]){
        System.out.println("enter first number of starting interval");
        double initial = sc.nextDouble();
        System.out.println("enter second number of the starting interval");
        double initial2 = sc.nextDouble();
        double width = initial2-initial;
        double mean=0,midvalue=0,sum=0,sumf=0;
        System.out.println("is the data continuous(Y:1,N:0)");
        switch(sc.nextInt()){
            case 1:
                midvalue = ((initial2+initial)/2);
                for(int i=0;i<values.length;i++){
                    sum+=values[i]*midvalue;
                    sumf+=values[i];
                    midvalue+=width;
                }
                mean=sum/sumf;
                break;
            case 0:
                System.out.println("enter first number of second interval");
                double initial3 = sc.nextDouble();
                double difference = (initial3-initial2)/2;
                initial-=difference;
                initial2+=difference;
                width = initial2-initial;
                midvalue = ((initial2+initial)/2);
                for(int i=0;i<values.length;i++){
                    sum+=values[i]*midvalue;
                    sumf+=values[i];
                    midvalue+=width;
                }
                mean=sum/sumf;
                break;
            default:
                System.out.println("invalid choice");
                System.exit(0);
        }
        return mean;
    }
    public static double[] median(double[] nums){
        double[] summary = new double[5];
        Arrays.sort(nums);
        int count = nums.length;
        summary[0] = nums[0];
        summary[4] = nums[count-1];
        if(Math.ceil(0.25*count)==(0.25*count)){
            summary[1]=(nums[(int)0.25*count]+nums[(int)(0.25*count)+1])/2;
        }else{
            summary[1]=nums[(int)Math.ceil(0.25*count)];
        }
        if(Math.ceil(0.5*count)==(0.5*count)){
            summary[2]=(nums[(int)0.5*count]+nums[(int)(0.5*count)+1])/2;
        }else{
            summary[2]=nums[(int)Math.ceil(0.5*count)];
        }
        if(Math.ceil(0.75*count)==(0.75*count)){
            summary[3]=(nums[(int)0.75*count]+nums[(int)(0.75*count)+1])/2;
        }else{
            summary[3]=nums[(int)Math.ceil(0.75*count)];
        }
        return summary;
    }
    public static double variances(double[] list,int c){
        double mean = nMean(list);
        double sum=0;
        for(int i=0;i<list.length;i++){
            sum+=Math.pow((list[i]-mean), 2);
        }
        if(c==0){
            return sum/(list.length-1);
        }else{
            return sum/(list.length);
        }
    }
    public static double standardDeviation(double[] list,int c){
        if(c==0){
            return Math.pow(variances(list, 0), 0.5);
        }else{
            return Math.pow(variances(list, 1), 0.5);  
        }
    }
    public static double covariances(double[] x,double[] y,int c){
        double xmean = nMean(x);
        double ymean = nMean(y);
        double sum=0;
        for(int i=0;i<x.length;i++){
            sum+=(xmean-x[i])*(ymean-y[i]);
        }
        if(c==0){
            return sum/(x.length-1);
        }else{
            return sum/(x.length);
        }
    }
    public static double correlation(double[] x,double[] y){
        double covariance = covariances(x, y, 0);
        double sx = standardDeviation(x, 0);
        double sy = standardDeviation(y, 0);
        return covariance/(sx*sy);
    }
    public static double pbscorrelation(double[] x0,double[] x1){
        double mean0 = nMean(x0);
        double mean1 = nMean(x1);
        double sx = standardDeviation(x0, 0);
        double p0 = x0.length/(x0.length+x1.length);
        double p1 = x1.length/(x0.length+x1.length);
        return ((mean0-mean1)/sx)*Math.pow((p0*p1), 0.5);
    }
    public static int permutations(int n,int r){
        return (factorial(n)/factorial(n-r));
    }
    public static int combinations(int n,int r){
        return (factorial(n)/(factorial(r)*factorial(n-r)));
    }
    public static double[] randomEV(int[] values,double[] probabilities){
        double expectation=0;
        for(int i=0;i<values.length;i++){
            expectation+=values[i]*probabilities[i];
        }
        double variance=0;
        for(int i=0;i<values.length;i++){
            variance+=Math.pow(values[i], 2)*probabilities[i];
        }
        variance-=Math.pow(expectation, 2);
        double[] out = {expectation,variance};
        return out;
    }
    public static double[] binomial(int trials,double probability){
        double[] distribution = new double[trials+1];
        for(int i=0;i<=trials;i++){
            distribution[i]=combinations(trials, i)*Math.pow(probability, i)*Math.pow(probability-1, trials-i);
        }
        return distribution;
    }
    public static double[] binomialEV(int trials,double probability){
        double[] out = {(double)trials*probability,(double)trials*probability*(1-probability)};
        return out;
    }
    public static double[] hypergeometric(int population,int sample,int trials){
        int max = Math.min(sample, population-sample);
        double[] distribution = new double[max+1];
        for(int i=0;i<=max;i++){
            distribution[i]=(combinations(sample, i)*combinations(population-sample, trials-i))/combinations(population, trials);
        }
        return distribution;
    }
    public static double[] hypergeometricEV(int population,int sample,int trials){
        double[] out = {(trials*sample)/population,((trials*sample)/population)*((population-sample)/population)*((population-trials)/(population-1))};
        return out;
    }
    public static double[] poisson(double rate){
        double[] distribution = new double[101];
        for(int i=0;i<101;i++){
            distribution[i]=(Math.pow(rate, i)*Math.pow(2.7182, (-1)*rate))/factorial(i);
        }
        return distribution;
    }
    public static int factorial(int n){
        int fact=1;
        for(int i=2;i<n;i++){
            fact*=i;
        }
        return fact;
    }
    public static void printlist(double[] a){
        System.out.print("[");
        for(int i=0;i<a.length;i++){
            System.out.print(a[i]+" ");
        }
        System.out.println("\b]");
    }
    public static double[] input(){
        ArrayList<Double> values = new ArrayList<Double>();
        Double n;
        while(true){
            n=sc.nextDouble();
            if(n==-1){
                break;
            }
            values.add(n);
        }
        double[] values1 = convert(values);
        return values1;
    }
    public static double[] convert(ArrayList<Double> list){
        double[] array = new double[list.size()];
        for(int i=0;i<list.size();i++){
            array[i]=list.get(i);
        }
        return array;
    }
}