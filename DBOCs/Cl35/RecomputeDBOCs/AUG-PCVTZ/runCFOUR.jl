point::Int64 = 2
# jobSplit::Int64 = 9
# numberOfJobs::Int64 = 1000

ROH::Float64=0.96000000                                                                  
ROCL::Float64=1.68300000                                                                 
A::Float64=107.00000

run(`xclean`)
# run(`./runGenerateCfourScript35.csh $(point) $(ROH) $(ROCL) $(A)`)
# run(`./runGenerateCfourScript37.csh $(point) $(ROH) $(ROCL) $(A)`)
# point += 1

rOH::Float64 = 0.6
while rOH < 2.8
    run(`xclean`)
    run(`./runGenerateCfourScript35.csh $(point) $(rOH) $(ROCL) $(A)`)
    run(`./runGenerateCfourScript37.csh $(point) $(rOH) $(ROCL) $(A)`)
    global rOH += 0.05
    global point += 1
end

rOCL::Float64 = 1.4
while rOCL < 3.5
    run(`xclean`)
    run(`./runGenerateCfourScript35.csh $(point) $(ROH) $(rOCL) $(A)`)
    run(`./runGenerateCfourScript37.csh $(point) $(ROH) $(rOCL) $(A)`)
    global rOCL += 0.05
    global point += 1
end

a::Float64 = 30
while rOCL < 180
    run(`xclean`)
    run(`./runGenerateCfourScript35.csh $(point) $(ROH) $(ROCL) $(a)`)
    run(`./runGenerateCfourScript37.csh $(point) $(ROH) $(ROCL) $(a)`)
    global a += 5
    global point += 1
end

# for i in 1:numberOfJobs
#     run(`qsub -e CH3OH_CBS_$(point)-$(point + jobSplit).e -o CH3OH_CBS_$(point)-$(point + jobSplit).o -l h_rt="11:59:00" -l mem=70G -l tmpfs=100G RunMolproJobs.csh $(point) $(point + jobSplit)`)
#     global point = point + jobSplit
#     println(point)
# end