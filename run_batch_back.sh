#!/bin/bash

steps=(     499    487    475    464    453    442    432    421    411	   401	  392	 382	373    365    355    347    338	   331	  323	 315	307    300    293    286    279	   272	  266	 259	253    247    241    235 )
steps=(     401    392   382	373    365    355    347    338	   331	  323	 315	307    300    293    286    279	   272	  266	 259	253    247    241    235 )
steps=(         235    247    241    253	   259	  266	 272	279    286    )

for step in "${steps[@]}";do
    echo "\n\n\n running ${step}..."
    ./galaticus_all_gal_pos.py params/batch/${step}.param
    ./main params/batch/${step}.param
done
