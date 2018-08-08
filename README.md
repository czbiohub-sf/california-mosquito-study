# skeeters

Files used to perform each analysis are stored in a separate sub-folder in the 'analysis' folder.

The 'data' folder contains small data files such as metadata information. Larger files such as sequence fastqs and fastas are stored on S3 (s3://czbiohub-mosquito/).

Figures generated from analyses are stored in the 'figures'directory.

The 'scripts' folder contains R, Python, Bash scripts etc. that could be re-used for different parts of the analyses. 

## How to contribute to this github repo

1. Open up a Terminal window on your local computer and go to the location where you want to save this repo. E.g.

```
cd /Users/lucy/projects
```

2. Clone the Github repo (i.e. download all the files and folders listed above to your local computer).

```
git clone https://github.com/czbiohub/skeeters.git
```

3. Now you should be able to see a 'skeeters' folder. To get the up-to-date files from this repo, 'pull' from this repo using

```
git pull
```

4. If you have edited the files on your local computer (or added/deleted files), you can upload the changes you have made to this Github repo. This involves first 'committing' the changes you have made, which means you select which of the changes you want to upload. And then 'pushing' those changes so that the Github repo contains your changes.

```
git commit . -m "commit message"
git push
```

