#!/bin/bash

steps=(     499    487    475    464    453    442    432    421    411	   401	  392	 382	373    365    355    347    338	   331	  323	 315	307    300    293    286    279	   272	  266	 259	253    247    241    235 )
stepsA=(    499  487  475  464  453  442  432  421  411	 401  392  )
stepsB=(    382  373  365  355  347  338  331  323  315  307  300  )
stepsC=(    293  286  279  272  266  259  253  247  241  235 )

for step in "${stepsB[@]}";do
    echo "\n\n\n running ${step}..."
    param=params/batch3/${step}.param
    ./galaticus_all_gal_pos_explicit_list.py ${param}
    ./main ${param}
done
