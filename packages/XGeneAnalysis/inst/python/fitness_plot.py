__author__ = 'DALE'

def ttest(a,b,c,d,e,f,equal_var=False):
    
    print('Debug')
    return ttest2(a,b,c,d,e,f,equal_var)

def fitness_plot(input_tsvfile='original_A.tsv',strain_name='A',gene_list=['YOR1','YCF1','SNQ2','PDR5','YBT1'],antibiotic_name='fluconazole',pvalue=0.05,dot_color_list=[(1,1,1)]*20,knockout_color=(0,0,0),line_color_list=[(0,0,1),(1,0.6,0),(0.5,0.5,0.5),(0.5,0.5,0.5)],*args,**kwargs):  #specific any gene knockout combinations with 16 antibiotic treatments (A,alpha)
    '''
    all the arguments:
    input_tsvfile, input file
    strain_name, this strain_name will be added to the plot title and file name
    gene_list, put all the genes in a list
    antibiotic_name, put the antibiotics of interest
    pvalue, set the threshold for t-test
    dot_color, input the colors of all the genes (wild type) in a list
    knockout_color, default black
    line_color_list, input three different colors in a list for increased, decreased, orphan, and nonsignificant lines
    *args and **kwargs, receive more incorrect arguments
    '''
    #process the input file and add a column of knouckout number
    from scipy.stats import ttest_ind_from_stats as ttest2
    import csv
    import matplotlib.pyplot as plt
    import numpy as np
    import os
    import pandas as pd
    

    file1=open(input_tsvfile,'r')
    file2=open('sort_info_'+strain_name+'_titled','w')
    reader=csv.reader(file1,delimiter='\t')
    writer=csv.writer(file2,delimiter='\t')
    for row in reader:
        if row[0]=='ID':
            writer.writerow(row[:18]+['Knockout_number']+row[18:])
        else:
            numknock=np.sum([int(i) for i in row[2:18]])
            writer.writerow(row[:18]+[numknock]+row[18:])
    file1.close()
    file2.close()

    #default parameters for figure
    gene_number=len(gene_list)
    side_number=gene_number/2
    if gene_number<4:
        figwidth=7.5
        space_len=0.1
        marker_size=7.75
        line_width=1.25
    elif gene_number>=4 and gene_number<=12:
        figwidth=1.25*gene_number+2.5
        if gene_number>10:
            marker_size=2.5
            line_width=0.5
            space_len=0.04
        else:
            marker_size=-0.875*gene_number+11.25
            line_width=-0.125*gene_number+1.75
            space_len=-0.01*gene_number+0.14
    elif gene_number>12:
        figwidth=17.5
        space_len=0.02
        marker_size=2.5
        line_width=0.5

    if gene_number<4:
        space_len=0.1
    elif gene_number>=4 and gene_number<=10:
        space_len=-0.01*gene_number+0.14
    elif gene_number>10:
        space_len=0.003*gene_number-0.007

    #create a pivot table by selected genes and antibiotics,get the mean,sum,std from the same genotype,output a tsv file.
    plt.figure(figsize=(figwidth,6))
    filetemp=open('sort_info_'+strain_name+'_mean','w')
    dataframe=pd.read_table('sort_info_'+strain_name+'_titled')
    table=pd.pivot_table(dataframe,values=antibiotic_name,index=gene_list,aggfunc=[np.mean,np.sum,np.std])
    table.to_csv('sort_info_'+strain_name+'_mean',sep='\t')

    #add a column of knockout number to the pivot table.
    fitness_list=[]
    file = open('sort_info_'+strain_name+'_mean_plot','w')
    writer= csv.writer(file,delimiter='\t')
    with open('sort_info_'+strain_name+'_mean','rb') as tsvfile:
        reader = csv.reader(tsvfile,delimiter='\t')
        for row in reader:
            if row[-1]=='std':
                continue
            if row[gene_number+2]=='':
                row[gene_number+2]=0
            knockout_number=0
            for i in range(gene_number):
               knockout_number = knockout_number + int(row[i])
            row.append(knockout_number)
            writer.writerow(row)
            fitness_list.append(float(row[gene_number]))
    file.close()
    tsvfile.close()

    #plot dots of even or odd number
    file = open('sort_info_'+strain_name+'_mean_plot_temp','w')
    writer= csv.writer(file,delimiter='\t')
    with open('sort_info_'+strain_name+'_mean_plot','rb') as tsvfile:
        reader = csv.reader(tsvfile,delimiter='\t')
        if gene_number%2==0:
            for row in reader:
                for i in range(gene_number):
                    if int(row[i])==0:
                        plt.plot(int(row[gene_number+3])-space_len*(side_number-1)+space_len*i,float(row[gene_number]),color=dot_color_list[i],marker='o',ms=marker_size)
                    else:
                        plt.plot(int(row[gene_number+3])-space_len*(side_number-1)+space_len*i,float(row[gene_number]),color=knockout_color,marker='o',ms=marker_size)
                writer.writerow(row)
        else:
            for row in reader:
                for i in range(gene_number):
                    if int(row[i])==0:
                        plt.plot(int(row[gene_number+3])-space_len*side_number+space_len*i,float(row[gene_number]),color=dot_color_list[i],marker='o',ms=marker_size)
                    else:
                        plt.plot(int(row[gene_number+3])-space_len*side_number+space_len*i,float(row[gene_number]),color=knockout_color,marker='o',ms=marker_size)
                writer.writerow(row)
    tsvfile.close()
    file.close()

    #plot lines and do t-test:1) one sample, grey line;2)no significance, grey line;3)increase;4)decrease
    if gene_number%2==0:
        with open('sort_info_'+strain_name+'_mean_plot','rb') as tsvfile:
            reader = csv.reader(tsvfile,delimiter='\t')
            for row1 in reader:
                file = open('sort_info_'+strain_name+'_mean_plot_temp','rb')
                reader2= csv.reader(file,delimiter='\t')
                for row2 in reader2:
                    sum=0
                    if int(row2[gene_number+3])-int(row1[gene_number+3])==1:
                        for i in range(gene_number):
                            sum+=abs(int(row1[i])-int(row2[i]))
                        if sum==1:
                            if float(row1[gene_number+2])==0 or float(row2[gene_number+2])==0:
                                plt.plot([int(row1[gene_number+3])+space_len*side_number+0.75*space_len,int(row2[gene_number+3])-space_len*side_number+0.25*space_len],[float(row1[gene_number]),float(row2[gene_number])],color=line_color_list[3],linewidth=line_width,ls='--')
                            else:
                                t_value,pval=ttest(float(row1[gene_number]),float(row1[gene_number+2]),float(row1[gene_number+1])/float(row1[gene_number]),
                                                                        float(row2[gene_number]),float(row2[gene_number+2]),float(row2[gene_number+1])/float(row2[gene_number]),
                                                                        equal_var=False)
                                #print t_value,pval
                                if pval>pvalue:
                                    plt.plot([int(row1[gene_number+3])+space_len*side_number+0.75*space_len,int(row2[gene_number+3])-space_len*side_number+0.25*space_len],[float(row1[gene_number]),float(row2[gene_number])],color=line_color_list[2],linewidth=line_width,ls='--')
                                else:
                                    if row2[gene_number]>row1[gene_number]:
                                        plt.plot([int(row1[gene_number+3])+space_len*side_number+0.75*space_len,int(row2[gene_number+3])-space_len*side_number+0.25*space_len],[float(row1[gene_number]),float(row2[gene_number])],color=line_color_list[0],linewidth=line_width)
                                    else:
                                        plt.plot([int(row1[gene_number+3])+space_len*side_number+0.75*space_len,int(row2[gene_number+3])-space_len*side_number+0.25*space_len],[float(row1[gene_number]),float(row2[gene_number])],color=line_color_list[1],linewidth=line_width)

    else:
        with open('sort_info_'+strain_name+'_mean_plot','rb') as tsvfile:
            reader = csv.reader(tsvfile,delimiter='\t')
            for row1 in reader:
                file = open('sort_info_'+strain_name+'_mean_plot_temp','rb')
                reader2= csv.reader(file,delimiter='\t')
                for row2 in reader2:
                    sum=0
                    if int(row2[gene_number+3])-int(row1[gene_number+3])==1:
                        for i in range(gene_number):
                            sum+=abs(int(row1[i])-int(row2[i]))
                        if sum==1:

                            if float(row1[gene_number+2])==0 or float(row2[gene_number+2])==0:
                                plt.plot([int(row1[gene_number+3])+space_len*side_number+0.75*space_len,int(row2[gene_number+3])-space_len*side_number-0.75*space_len],[float(row1[gene_number]),float(row2[gene_number])],color=line_color_list[3],linewidth=line_width,ls='--')
                            else:
                                t_value,pval=ttest(float(row1[gene_number]),float(row1[gene_number+2]),float(row1[gene_number+1])/float(row1[gene_number]),
                                                                        float(row2[gene_number]),float(row2[gene_number+2]),float(row2[gene_number+1])/float(row2[gene_number]),
                                                                        equal_var=False)

                                if pval>pvalue:
                                    plt.plot([int(row1[gene_number+3])+space_len*side_number+0.75*space_len,int(row2[gene_number+3])-space_len*side_number-0.75*space_len],[float(row1[gene_number]),float(row2[gene_number])],color=line_color_list[2],linewidth=line_width,ls='--')
                                else:
                                    if row2[gene_number]>row1[gene_number]:
                                        plt.plot([int(row1[gene_number+3])+space_len*side_number+0.75*space_len,int(row2[gene_number+3])-space_len*side_number-0.75*space_len],[float(row1[gene_number]),float(row2[gene_number])],color=line_color_list[0],linewidth=line_width)
                                    else:
                                        plt.plot([int(row1[gene_number+3])+space_len*side_number+0.75*space_len,int(row2[gene_number+3])-space_len*side_number-0.75*space_len],[float(row1[gene_number]),float(row2[gene_number])],color=line_color_list[1],linewidth=line_width)

    #plot x,y axis label,title,set limits for x,y axis
    plt.ylabel('fitness value')
    plt.xlabel('knockout number')
    plt.title('Fitness landscape of '+strain_name+' strain for '+antibiotic_name)
    plt.xlim(-0.6,gene_number+0.6)
    plt.ylim(min(fitness_list)-(max(fitness_list)-min(fitness_list))*0.05,(max(fitness_list)-min(fitness_list))*0.05+max(fitness_list))

    file.close()
    tsvfile.close()

    #plot the legend, use the wild type fitness to determine where the legend should be placed
    with open('sort_info_'+strain_name+'_mean_plot_temp','r') as file:
        reader=csv.reader(file,delimiter='\t')
        for row in reader:
            if row[:gene_number]==['0']*gene_number:
                value_range= abs(min(fitness_list)-(max(fitness_list)-min(fitness_list))*0.05-(max(fitness_list)-min(fitness_list))*0.05-max(fitness_list))
                value_to_top=(max(fitness_list)-min(fitness_list))*0.05+max(fitness_list)-float(row[gene_number])
                percentage=value_to_top/value_range
                #print value_range,value_to_top,percentage
                #Checks to see how much of the top left of the plot is free
                if percentage>0.22:
                    for i in range(gene_number):
                        plt.plot(-0.4+0.1*i,max(fitness_list)-0.01*(max(fitness_list)-min(fitness_list)),color=dot_color_list[i],marker='o',markersize=6)
                        plt.text(-0.45+0.1*i,max(fitness_list)-0.04*(max(fitness_list)-min(fitness_list)),gene_list[i],rotation=-50,size=6.5)

                    plt.plot(-0.4+0.1*(i+1),max(fitness_list)-0.01*(max(fitness_list)-min(fitness_list)),'ko')
                    plt.text(-0.45+0.1*(i+1),max(fitness_list)-0.04*(max(fitness_list)-min(fitness_list)),'Knockout',rotation=-50,size=6)

                    line_list=['increase','decrease','non.sig']
                    for i in range(3):
                        plt.plot([-0.4,-0.2],[max(fitness_list)-(0.11+0.03*i)*(max(fitness_list)-min(fitness_list)),max(fitness_list)-(0.11+0.03*i)*(max(fitness_list)-min(fitness_list))],color=line_color_list[i])
                        plt.text(-0.15,max(fitness_list)-(0.115+0.03*i)*(max(fitness_list)-min(fitness_list)),line_list[i],fontsize=7.5)
                else:
                    for i in range(gene_number):
                        plt.plot(-0.4+0.1*i,max(fitness_list)-0.8*(max(fitness_list)-min(fitness_list)),color=dot_color_list[i],marker='o',markersize=6)
                        plt.text(-0.45+0.1*i,max(fitness_list)-0.83*(max(fitness_list)-min(fitness_list)),gene_list[i],rotation=-50,size=6.5)

                    plt.plot(-0.4+0.1*(i+1),max(fitness_list)-0.8*(max(fitness_list)-min(fitness_list)),'ko')
                    plt.text(-0.45+0.1*(i+1),max(fitness_list)-0.83*(max(fitness_list)-min(fitness_list)),'Knockout',rotation=-50,size=6)

                    line_list=['increase','decrease','non.sig']
                    for i in range(3):
                        plt.plot([-0.4,-0.2],[max(fitness_list)-(0.89+0.03*i)*(max(fitness_list)-min(fitness_list)),max(fitness_list)-(0.89+0.03*i)*(max(fitness_list)-min(fitness_list))],color=line_color_list[i])
                        plt.text(-0.15,max(fitness_list)-(0.895+0.03*i)*(max(fitness_list)-min(fitness_list)),line_list[i],fontsize=7.5)

    #remove the temporary files
    os.remove('sort_info_'+strain_name+'_mean_plot_temp')
    os.remove('sort_info_'+strain_name+'_mean')
    os.remove('sort_info_'+strain_name+'_mean_plot')
    os.remove('sort_info_'+strain_name+'_titled')
    
    #the format of file saving
    plt.savefig(strain_name+'_%dgenes_'%gene_number+antibiotic_name+'.pdf',format='pdf',dpi=900)
