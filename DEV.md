## development

grep -r "src.main.python" . | grep ".py" | grep -v "pycache" | cut -d ":" -f 1 | uniq | grep polygenic | xargs -i bash -c "sed -i 's/src.main.python/polygenic/g' {}"

### test
```
python3 setup.py nosetests -s
pytest -s -k "testPolygenicForGbeModel"
```
### build
```
python3 setup.py sdist bdist_wheel
```
### upload
```
twine upload --repository testpypi dist/polygenic-0.0.1*
twine upload --repository pypi dist/polygenic-1.0.0*
```
### install
```
pip3 install --upgrade --no-deps polygenic==1.1.5
```
