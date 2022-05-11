cd ~/dct-policy
source venv/bin/activate
# I had to install python3.6 separately, in case
# you have this version installed you may use
# ```python run.py < $1```
# instead
python3.6 run_col.py < $1
