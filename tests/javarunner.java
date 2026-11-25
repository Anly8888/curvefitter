import ij.measure.CurveFitter;
import ij.measure.Minimizer;
import java.util.*;

// Compiled against the pinned, unmodified ImageJ sources; emits committed oracle fixtures.
class JavaRunner {
    static String number(double d) { return Double.isFinite(d)?Double.toString(d):"null"; }
    static String array(double[] a) {
        if (a==null) return "null";
        StringJoiner out=new StringJoiner(",","[","]");
        for(double d:a)out.add(number(d));
        return out.toString();
    }
    static void emitFit(String name,int type,double[] x,double[] y,double[] initial,double[] weights) {
        CurveFitter f=new CurveFitter(x,y);
        f.getMinimizer().setMaximumThreads(1);
        f.getMinimizer().setRandomSeed(0);
        if(initial!=null)f.setInitialParameters(initial.clone());
        if(weights!=null)f.setWeights(weights);
        f.doFit(type);
        System.out.println("{\"name\":\""+name+"\",\"type\":"+type+",\"x\":"+array(x)+",\"y\":"+array(y)
            +",\"initial\":"+array(initial)+",\"weights\":"+array(weights)+",\"params\":"+array(f.getParams())
            +",\"status\":"+f.getStatus()+",\"iterations\":"+f.getIterations()+",\"sse\":"+number(f.getSumResidualsSqr())
            +",\"r2\":"+number(f.getRSquared())+",\"sd\":"+number(f.getSD())+",\"goodness\":"+number(f.getFitGoodness())
            +",\"residuals\":"+array(f.getResiduals())+"}");
    }
    static double[] params(int type) {
        switch(type){
            case 0:return new double[]{2,3};case 1:return new double[]{2,3,.2};case 2:return new double[]{2,3,.2,-.1};case 3:return new double[]{2,3,.2,-.1,.05};
            case 4:case 15:return new double[]{2,.3};case 5:case 16:return new double[]{2,1.3};case 6:return new double[]{2,3};
            case 7:return new double[]{8,2,3,1};case 8:return new double[]{.1,2,2,1};case 9:return new double[]{2,3,-1};
            case 10:case 14:return new double[]{0,2,3,8};case 11:case 13:return new double[]{2,.7,1};case 12:return new double[]{1,4,2,.7};
            case 21:return new double[]{3,2,.7};case 22:return new double[]{3,.7};case 23:return new double[]{3,.7,1.5};case 24:return new double[]{2,3,2,.8};
            default:double[] p=new double[CurveFitter.getNumParams(type)];Arrays.fill(p,.01);p[0]=2;p[1]=1;return p;
        }
    }
    static void emitMin(String name,int mode) {
        Minimizer m=new Minimizer();m.setMaximumThreads(1);m.setRandomSeed(31);
        final int[] calls={0};
        m.setFunction((p,x)->{
            if(mode==4 && ++calls[0]>100)return Double.NaN;
            if(mode==1)return Double.NaN;
            if(mode==2){m.abort();return p[0]*p[0];}
            return (p[0]-2)*(p[0]-2)+(p[1]+3)*(p[1]+3)+(mode==5?m.getCompletedMinimizations():0);
        },2);
        if(mode==3)m.setMaxIterations(1);
        if(mode==5)m.setMaxRestarts(1);
        m.setMaxError(1e-10,1e-12);
        int status=mode==4?m.minimizeOnce(new double[]{0,0},new double[]{1,1}):m.minimize(new double[]{0,0},new double[]{1,1});
        System.out.println("{\"name\":\""+name+"\",\"minimizer\":"+mode+",\"params\":"+array(m.getParams())
            +",\"status\":"+status+",\"iterations\":"+m.getIterations()+",\"completed\":"+m.getCompletedMinimizations()+"}");
    }
    public static void main(String[] args) {
        for(int type=0;type<25;++type){
            double[] x=new double[40],y=new double[40],p=params(type),w=new double[40];
            for(int i=0;i<x.length;++i){x[i]=.2+4.*i/39;y[i]=CurveFitter.f(type,p,x[i])+.001*Math.sin(i*1.7);w[i]=i%5==0?0:1+i%3;}
            emitFit("default_"+type,type,x,y,null,null);
            emitFit("seeded_"+type,type,x,y,p,null);
            emitFit("weighted_"+type,type,x,y,p,w);
        }
        emitFit("constant_x",0,new double[]{1,1,1},new double[]{2,3,4},null,null);
        emitFit("negative_exp",15,new double[]{1,2,3},new double[]{-2,-4,-8},null,null);
        emitFit("power_origin",16,new double[]{0,1,2,3},new double[]{0,2,8,18},null,null);
        emitFit("nan_filter",0,new double[]{0,Double.NaN,1},new double[]{1,99,3},null,null);
        emitFit("power_origin_negative",16,new double[]{0,1,2,3},new double[]{0,2,1,.5},null,null);
        for(int mode=0;mode<6;++mode)emitMin("minimizer_"+mode,mode);
    }
}
