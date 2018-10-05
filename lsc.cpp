#include<iostream>
#include<fstream>
#include<string>
#include<istream>

using namespace std;
 
main(){

ifstream datain;
ofstream dataout;

string inputname;
string outputname;
int count=0;
double m1,m2;
double *col1, *col2;
string line; 
double temp;

cout<<"This program is special for *n by 2* datafile"<<endl;
cout<<"Please indicate the file to be read: "<<endl;
cin>> inputname;
datain.open(inputname.c_str());

while(getline(datain, line))
	count ++;
datain.close();

count=(count-4)/3*2;

datain.open(inputname.c_str());

col1 = new double[count]();
col2 = new double[count]();



/*cout <<"There are *"<<count<<"* lines in the file"<<endl;

cout << "(data_input) * (multiplier) = data_output"<<endl;
cout << "Enter the multiplier for the FIRST column"<< endl;
cin>> m1;
cout << "Enter the multiplier for the SECOND column"<< endl;
cin>> m2; */

cout<<"Please name the outputfile: "<<endl;
cin >> outputname;
dataout.open(outputname.c_str());

datain.ignore(1000, '\n');
datain.ignore(1000, '\n');
datain.ignore(1000, '\n');
datain.ignore(1000, '\n');


for(int i=0;i<count;i++)
 {
   getline(datain ,line);
   dataout<<line<<endl;  

 }
/*for(int i=0;i<count;i++){
  datain >> col1[i];
  datain >> col2[i];
}

for(int i=0;i<count;i++){
  dataout<< col1[i];
  dataout<< "    ";
  dataout<< col2[i];
  dataout<<endl;
}
*/

delete[] col1;
delete[] col2;


datain.close();
dataout.close();
}
