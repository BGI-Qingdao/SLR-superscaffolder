# How to run this test case . 

## step 0 , install software

```
YOUR-DOWNLOAD-DIR/install.sh YOUR-INSTALL-DIR
```

## step 1 , edit conf.ini

There are 2 variabies that need to be assigned :

```
STLFR_ASSEMBLER_DIR=YOUR-INSTALL-DIR         # stLFR Scaffold Assembler installation directory
BWA_DIR=YOUR-BWA-DIR
```

## step 2 , prepare scripts

```
YOUR-INSTALL-DIR/scaffold/prepare.sh conf.ini
```
Above command will create a new folder *scaff_example* and put all scripts into it.

## step 3 , run scaffold scripts

```
cd scaff_example && ./run.sh
```

## done !
