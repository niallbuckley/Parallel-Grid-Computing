while read -u 3 name; do
  ssh -o StrictHostKeyChecking=no ${name%.ucc.ie} "echo done $name"
done 3< machines
