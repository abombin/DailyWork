aws s3api list-buckets --query "Buckets[].Name" >bucket.list

for i in $(cat bacteria.list)
do
	if grep -qw "$i"  bucket.list; then
		aws s3 cp ./bactopia_gtdbtk/"$i"/"$i"*.zip s3://"$i"
	else
		aws s3 mb s3://"$i"
		aws s3 cp ./bactopia_gtdbtk/"$i"/"$i"*.zip s3://"$i"
	fi
done

# transfer metadata to dynamo db
aws s3 rm s3://bacteria-metadata/metaDat.csv
sleep 60s
aws s3 cp ./metaData/metaDat.csv s3://bacteria-metadata
