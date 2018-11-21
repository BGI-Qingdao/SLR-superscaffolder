# stLFR denovo prototype project

## Install
    
    ./intall.sh dir_to_install
    

## pipeline 


### common

####	拷贝对应流程的conf.ini 文件到项目目录下,并根据自己的配置修改此文件

>cp /home/software/stLFR/script/MST/conf.ini /home/project/xxx/

####	生成流程

> cd  /home/project/xxx/
 
> /home/software/stLFR/script/MST/prepare.sh ./conf.ini

#### 运行脚本。 上一步根据配置，会生成项目目录/home/project/xxx/test， 内部有全部脚本。

> cd test

> ./run.sh # 或者按需要单步执行具体脚本。


### part 1 . contig extern by PE and barcode info and DBG .

TODO

### part 2 . contig extern by barcode and DBG .

TODO 

### part 3 . scaffolding by barcode .

TODO 

