#!/usr/bin/zsh

####################
## User parameter ##
###################################

user=HDB
project=TnpBPE
pam_type=Cpf1
pam_pos=Forward
thread=16

gap_open=-10 ## default
gap_extend=1 ## default

###################################
# TODO: if there is space characters...
while read python_path; do
    python=$python_path
done <../PythonPath.txt
done <../PythonPath.txt

[ ! -d ./Output/${user} ] && eval mkdir ./Output/${user}
[ ! -d ./Output/${user}/${project} ] && eval mkdir ./Output/${user}/${project}
[ ! -d ./Output/${user}/${project}/Log ] && eval mkdir ./Output/${user}/${project}/Log

nohup $python ./Run_indel_searcher.py --python $python --user $user --project $project --pam_type $pam_type --pam_pos $pam_pos -t $thread >./Output/${user}/${project}/Log/log.txt 2>&1
