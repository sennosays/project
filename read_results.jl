

fp = open(name,"r"); 
this_length = read(fp,Int,1)[1]; 
density = read(fp,Float64,this_length);
close(fp); 

