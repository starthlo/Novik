# Generated migration for profile_completed field

from django.db import migrations, models


class Migration(migrations.Migration):
    dependencies = [
        ("api", "0006_conversation_delete_patientcontext"),
    ]

    operations = [
        migrations.AddField(
            model_name="customuser",
            name="profile_completed",
            field=models.BooleanField(
                default=False,
                help_text="Indicates if the user has completed their profile information",
            ),
        ),
    ]
