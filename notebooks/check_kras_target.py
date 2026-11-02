from chembl_webresource_client.new_client import new_client

target = new_client.target
res = target.filter(pref_name__icontains='KRAS')
for t in res:
    print(t['target_chembl_id'], '|', t['pref_name'], '|', t['organism'], '|', t['target_type'])

print('---')
res2 = target.filter(pref_name__icontains='KRas')
for t in res2:
    print(t['target_chembl_id'], '|', t['pref_name'], '|', t['organism'], '|', t['target_type'])
