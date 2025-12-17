#include<bits/stdc++.h>
using namespace std;
 
void printfunc(vector<int>cof)
{
    //reverse(cof.begin(),cof.end());
    int n=cof.size();
    for(int i=0; i<cof.size(); i++)
    {
        if(cof[i]>0&&i!=0)cout<<'+';
        if(cof[i]!=0)
        {
            cout<<cof[i];
            if(n-i-1>1)
            {
                cout<<"x^"<<n-i-1;
            }
            else if(i==n-1)
            {
                //cout<<cof[i]<<endl;
            }
            else
            {
                cout<<"x";
                //cout<<endl;
            }
        }
 
 
    }
    cout<<endl;
}
 
float func(vector<int>cof,float x)
{
    float ret=0;
    float xp=1;
    reverse(cof.begin(),cof.end());
 
    for(int i=0; i<cof.size(); i++)
    {
        float temp=xp*cof[i];
        ret=ret+temp;
        xp=xp*x;
    }
    return ret;
}
int main()
{
    //cout<<"Hello world\n";
    int n;
    cin>>n;
    vector<int>cof;
    for(int i=0; i<=n; i++)
    {
        int x;
        cin>>x;
        cof.push_back(x);
    }
    printfunc(cof);
 
    float itstart=1000000.355;
 
    vector<float>vv;
   /* for(float i=-itstart; i<itstart; i+=.5)
    {
        float it2=i+.5;
        float f1=func(cof,i);
        float f2=func(cof,it2);
        if(f1*f2<0){
            vv.push_back(i);
        }
    }*/
    if(vv.size()==0){
        vv.push_back(0);
    }
    else cout<<vv.size()<<endl;
 
    vector<float>res;
    float E=.00001;
 
    for(int i=0;i<vv.size();i++){
            float it=vv[i];
            float it2=it+.5;
 
        cout<<"Iteration starts for :"<<it<< " "<<it2<<endl;
 
        int cnt=1;
        map<float,bool>mp;
 
        while(1){
                cout<<"Number of iteration :"<<cnt<<endl;
 
            //1,2
            float f1=func(cof,it);
            float f2=func(cof,it2);
            // calculate 3
            float xr=it2-f2*(it-it2)/(f1-f2);
            cout<<it<<" "<<it2<<" "<<xr<<endl;
 
            float tmp=xr-it2;
            if(tmp<0)tmp=-tmp;
            tmp=tmp/xr;
 
 
            if(((tmp-E)<0)||mp[xr]){
                    res.push_back(xr);
                    break;
 
            }
            it=it2;
            it2=xr;
            mp[xr]=1;
            cnt++;
 
        }
    }
 
    cout<<"Results are : "<<endl;
    for(auto i:res){
        cout<<i<<" ";
    }
    cout<<endl;
 
}